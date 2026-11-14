#!/usr/bin/env python3
"""
ContextClimate hourly archive builder.

Builds one compact hourly record per station in the shared roster
(public/station-roster.js in josht3113/context-climate) from IEM's ASOS
archive, for the Streak Tracker. The tool reads these files instead of
pulling decades of observations from IEM on every visit, and only asks IEM
for the years after the archive's `through` year.

Output (under --out):
  hourly/<ICAO>/<decade>.bin.gz   one file per station per decade
  index.json                      stations, decades, `through`, variable list

Hourly values follow the Streak Tracker's rules exactly (see parse_station):
times are local standard time (LST), and each LST hour is one cell.

  tmpf dwpf feel   state   hundredths of a degree F
  vsby             state   hundredths of a mile
  sky              state   oktas x2 (highest layer: CLR 0, FEW 3, SCT 7, BKN 13, OVC/VV 16)
  wind             state   mph, rounded
  gust             event   mph; highest gust or PK WND in the clock hour; 0 = wind reported, no gust
  p01i             event   hundredths of an inch from the routine ob; -1 = trace
  wx               event   weather bits of every ob in the hour (rain 1, snow 2, fzra 4,
                           sleet 8, thunder 16, fog 32, smoke 64, haze 128)

State cells cover (H-1):45 to H:44; the routine ob wins and a special fills a
cell only when the routine value is missing. Missing = -32768.

File layout (little endian, then gzip):
  b'CCHR' | u8 version=2 | u8 nvars | u16 first_year | u16 n_years
  | nvars x (u8 variable index into VARS, u8 transform)
  | per variable: n_hours x int16, byte-shuffled (all low bytes, then all high bytes)
  n_hours = sum of the LST years' lengths (8760 / 8784), years consecutive.
  Transforms (lossless; missing stays -32768):
    0 raw
    1 delta from the previous non-missing value (first from 0)
    2 difference from tmpf in the same hour (only used when tmpf is present
      wherever the variable is; otherwise 1)

Modes:
  full         every decade for every station (cold start)
  incremental  the decade containing `through` for every station, plus a full
               history for any roster station not yet in the index
  stations     every decade for --stations
  index        rebuild index.json from the files under --out (no fetching)
"""
import argparse, datetime as dt, gzip, http.client, json, math, os, re, socket, struct, sys, time
import urllib.error, urllib.request
from zoneinfo import ZoneInfo

VARS = ['tmpf', 'dwpf', 'feel', 'vsby', 'sky', 'wind', 'gust', 'p01i', 'wx']
STATE_VARS = ['tmpf', 'dwpf', 'feel', 'vsby', 'sky', 'wind']
MISSING = -32768
TRACE_CODE = -1
FIELDS = ['tmpf', 'dwpf', 'feel', 'vsby', 'skyc1', 'skyc2', 'skyc3', 'skyc4', 'sknt', 'gust',
          'peak_wind_gust', 'peak_wind_time', 'p01i', 'metar']   # metar must stay last
HOUR = 3600000
KT_TO_MPH = 1.15078
ROUTINE_MIN_OBS = 200
YEAR_FLOOR = 1928
STATIONS_PER_REQUEST = 10
REQUEST_GAP_S = 1.5
RETRY_WAITS = (10, 30, 90)
ROSTER_URL = 'https://raw.githubusercontent.com/josht3113/context-climate/main/public/station-roster.js'
IEM_URL = 'https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py'
WX_BITS = {'rain': 1, 'snow': 2, 'fzra': 4, 'sleet': 8, 'ts': 16, 'fog': 32, 'smoke': 64, 'haze': 128}


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


# ── JavaScript-compatible helpers ───────────────────────────────────────────
def js_round(x):                       # Math.round: halves go up
    return math.floor(x + 0.5)


_NUM = re.compile(r'^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?')
def num_or_none(s):                    # numOrNull: trim, '' / 'M' -> None, parseFloat
    if s is None:
        return None
    s = s.strip()
    if not s or s == 'M':
        return None
    m = _NUM.match(s)
    return float(m.group(0)) if m else None


_off_cache = {}
def std_offset_min(iana, year):
    """Standard offset = the smallest offset seen across the year (1974 kept DST all winter)."""
    k = (iana, year)
    if k not in _off_cache:
        try:
            z = ZoneInfo(iana)
            v = min(int(dt.datetime(year, m, 15, 12, tzinfo=dt.timezone.utc).astimezone(z).utcoffset().total_seconds() // 60)
                    for m in range(1, 13))
        except Exception:
            v = -300
        _off_cache[k] = v
    return _off_cache[k]


def year_of_ms(ms):
    return (dt.datetime(1970, 1, 1) + dt.timedelta(milliseconds=ms)).year


def year_start_ms(y):
    return int((dt.datetime(y, 1, 1) - dt.datetime(1970, 1, 1)).total_seconds() * 1000)


def to_lst(s, iana):
    if not s or len(s) < 16:
        return None
    try:
        uy, um, ud, uh, umi = int(s[0:4]), int(s[5:7]), int(s[8:10]), int(s[11:13]), int(s[14:16])
        t = dt.datetime(uy, um, ud, uh, umi)
    except ValueError:
        return None
    return int((t - dt.datetime(1970, 1, 1)).total_seconds() * 1000) + std_offset_min(iana, uy) * 60000


def sky_code_to_oktas2(code):
    return {'CLR': 0, 'SKC': 0, 'NSC': 0, 'NCD': 0, 'FEW': 3, 'SCT': 7, 'BKN': 13, 'OVC': 16, 'VV': 16}.get(code)


_WX_TOKEN = re.compile(r'^(\+|-|VC)?(MI|PR|BC|DR|BL|SH|TS|FZ)?((?:DZ|RA|SN|SG|IC|PL|IP|GR|GS|UP|BR|FG|FU|VA|DU|SA|HZ|PY|PO|SQ|FC|SS|DS)*)$')
def wx_bits(metar):
    """The Streak Tracker's classifier: present-weather tokens in the body before RMK."""
    if not metar:
        return None
    up = metar.strip().upper()
    if not up or up == 'M':
        return None
    r = up.find(' RMK')
    if r >= 0:
        up = up[:r]
    bits = 0
    for t in up.split():
        if len(t) < 2:
            continue
        m = _WX_TOKEN.match(t)
        if not m:
            continue
        pre, desc, codes = m.group(1) or '', m.group(2) or '', m.group(3) or ''
        if not desc and not codes:
            continue
        cs = {codes[i:i + 2] for i in range(0, len(codes), 2)}
        liquid = 'RA' in cs or 'DZ' in cs
        if desc == 'TS': bits |= WX_BITS['ts']
        if desc == 'FZ' and liquid: bits |= WX_BITS['fzra']
        if desc != 'FZ' and pre != 'VC' and liquid: bits |= WX_BITS['rain']
        if 'SN' in cs or 'SG' in cs: bits |= WX_BITS['snow']
        if 'PL' in cs or 'IP' in cs: bits |= WX_BITS['sleet']
        if 'FG' in cs: bits |= WX_BITS['fog']
        if 'FU' in cs: bits |= WX_BITS['smoke']
        if 'HZ' in cs: bits |= WX_BITS['haze']
    return bits


def mode_of_hist(h):
    best, best_n = -1, 0
    for m in range(60):
        if h[m] > best_n:
            best_n, best = h[m], m
    return best


# ── Hourly selection (port of the Streak Tracker's parseChunk, all groups at once) ──
CLAMPED = {}
def _cl(v, name):
    if v is None or -32767 <= v <= 32767:
        return v
    CLAMPED[name] = CLAMPED.get(name, 0) + 1
    return 32767 if v > 0 else -32767


def parse_station(lines, idx, iana, y1, y2):
    """lines: CSV rows (already split) for ONE station, time order. idx: column index map.
    Returns {Y: None | {var: list[int]}} for LST years y1..y2 (None = no obs at all)."""
    rows = []
    hist_by_year, hist_all = {}, [0] * 60
    pcp_seen = set()
    sky_cols = [idx[f] for f in ('skyc1', 'skyc2', 'skyc3', 'skyc4') if f in idx]
    im = idx['metar']
    for p in lines:
        lst = to_lst(p[idx['valid']], iana)
        if lst is None:
            continue
        d = dt.datetime(1970, 1, 1) + dt.timedelta(milliseconds=lst)
        mi, ly = d.minute, d.year
        h = hist_by_year.setdefault(ly, [0] * 60)
        h[mi] += 1; hist_all[mi] += 1
        r = {'lst': lst, 'mi': mi, 'ly': ly}
        g = lambda f: p[idx[f]] if idx.get(f, -1) >= 0 and idx[f] < len(p) else ''
        for f in ('tmpf', 'dwpf', 'feel', 'vsby'):
            v = num_or_none(g(f))
            r[f] = None if v is None else _cl(js_round(v * 100), f)
        sk = None
        for c in sky_cols:
            code = (p[c] if c < len(p) else '').strip()
            if not code or code == 'M':
                continue
            o = sky_code_to_oktas2(code)
            if o is not None and (sk is None or o > sk):
                sk = o
        r['sky'] = sk
        k = num_or_none(g('sknt'))
        r['wind'] = None if k is None else _cl(js_round(k * KT_TO_MPH), 'wind')
        gu = num_or_none(g('gust'))
        r['seen'] = r['wind'] is not None or gu is not None
        r['gust'] = None if gu is None else _cl(js_round(gu * KT_TO_MPH), 'gust')
        pk = num_or_none(g('peak_wind_gust'))
        if pk is not None:
            pt = to_lst(g('peak_wind_time').strip(), iana) if 'peak_wind_time' in idx else None
            if not (pt is not None and pt <= lst and lst - pt <= 2 * HOUR):
                pt = lst - 30 * 60000
            r['pk'] = _cl(js_round(pk * KT_TO_MPH), 'gust'); r['pk_lst'] = pt
        r['bits'] = wx_bits(','.join(p[im:]))
        s = g('p01i').strip()
        if s == 'T':
            r['p'] = TRACE_CODE
        else:
            v = num_or_none(s)
            r['p'] = None if v is None else _cl(js_round(v * 100), 'p01i')
        if r['p'] is not None:
            pcp_seen.add(year_of_ms((lst // HOUR) * HOUR))
        rows.append(r)

    global_r = mode_of_hist(hist_all)
    routine_of = {}
    for yr, h in hist_by_year.items():
        routine_of[yr] = mode_of_hist(h) if sum(h) >= ROUTINE_MIN_OBS else global_r

    def dist(r):
        R = routine_of.get(r['ly'], global_r)
        dd = abs(r['mi'] - R)
        return min(dd, 60 - dd)

    bufs = {}
    def buf_for(Y):
        b = bufs.get(Y)
        if b is None:
            n = (year_start_ms(Y + 1) - year_start_ms(Y)) // HOUR
            b = {'vals': {v: [MISSING] * n for v in VARS}, 'dist': {v: [255] * n for v in STATE_VARS + ['p01i']}}
            bufs[Y] = b
        return b

    def clock_cell(ms):
        cell = (ms // HOUR) * HOUR
        Y = year_of_ms(cell)
        if Y < y1 or Y > y2:
            return None
        return buf_for(Y), (cell - year_start_ms(Y)) // HOUR

    for r in rows:
        d = dist(r)
        cell = ((r['lst'] + 15 * 60000) // HOUR) * HOUR
        Y = year_of_ms(cell)
        if y1 <= Y <= y2:
            b = buf_for(Y)
            i = (cell - year_start_ms(Y)) // HOUR
            for v in STATE_VARS:
                if r[v] is not None and d < b['dist'][v][i]:
                    b['dist'][v][i] = d; b['vals'][v][i] = r[v]
        # gust: highest in the clock hour (MISSING compares below every value)
        if r['seen']:
            c = clock_cell(r['lst'])
            if c:
                a = c[0]['vals']['gust']; val = 0 if r['gust'] is None else r['gust']
                if a[c[1]] < val: a[c[1]] = val
        if 'pk' in r:
            c = clock_cell(r['pk_lst'])
            if c:
                a = c[0]['vals']['gust']
                if a[c[1]] < r['pk']: a[c[1]] = r['pk']
        if r['bits'] is not None:
            c = clock_cell(r['lst'])
            if c:
                a = c[0]['vals']['wx']; prev = a[c[1]]
                a[c[1]] = (0 if prev == MISSING else prev) | r['bits']
        c = clock_cell(r['lst'])
        if c:
            b, i = c
            if not (d > 10 or d >= b['dist']['p01i'][i]):
                b['dist']['p01i'][i] = d
                b['vals']['p01i'][i] = 0 if r['p'] is None else r['p']

    out = {}
    for Y in range(y1, y2 + 1):
        b = bufs.get(Y)
        if b is None:
            out[Y] = None
            continue
        if Y not in pcp_seen:
            b['vals']['p01i'] = [MISSING] * len(b['vals']['p01i'])
        out[Y] = b['vals']
    return out


# ── Encoding ────────────────────────────────────────────────────────────────
def year_hours(y):
    return (year_start_ms(y + 1) - year_start_ms(y)) // HOUR


TRANSFORM = {'tmpf': 1, 'dwpf': 1, 'feel': 2}


def _delta(vals):
    out, prev = [], 0
    for v in vals:
        if v == MISSING:
            out.append(MISSING)
        else:
            out.append(v - prev); prev = v
    return out


def _undelta(vals):
    out, prev = [], 0
    for v in vals:
        if v == MISSING:
            out.append(MISSING)
        else:
            prev += v; out.append(prev)
    return out


def encode(first_year, years):
    """years: list (consecutive from first_year) of None | {var: list[int]}."""
    planes = {}
    for v in VARS:
        vals = []
        for k, yv in enumerate(years):
            n = year_hours(first_year + k)
            vals.extend([MISSING] * n if yv is None else yv[v])
        planes[v] = vals
    out = bytearray(b'CCHR')
    out += struct.pack('<BBHH', 2, len(VARS), first_year, len(years))
    body = bytearray()
    for vi, v in enumerate(VARS):
        vals, tr = planes[v], TRANSFORM.get(v, 0)
        if tr == 2 and any(x != MISSING and t == MISSING for x, t in zip(vals, planes['tmpf'])):
            tr = 1
        if tr == 1:
            vals = _delta(vals)
        elif tr == 2:
            vals = [MISSING if x == MISSING else x - t for x, t in zip(vals, planes['tmpf'])]
        if any(not (-32767 <= x <= 32767) for x in vals if x != MISSING):
            tr, vals = 0, planes[v]
        out += bytes([vi, tr])
        raw = struct.pack('<%dh' % len(vals), *vals)
        body += raw[0::2] + raw[1::2]
    return gzip.compress(bytes(out + body), compresslevel=9, mtime=0)


def decode(blob):
    if blob[:2] == b'\x1f\x8b':
        blob = gzip.decompress(blob)
    assert blob[:4] == b'CCHR', 'bad magic'
    ver, nv, fy, ny = struct.unpack('<BBHH', blob[4:10])
    assert ver == 2, f'unsupported version {ver}'
    spec = [(blob[10 + 2 * i], blob[11 + 2 * i]) for i in range(nv)]
    pos = 10 + 2 * nv
    n = sum(year_hours(fy + k) for k in range(ny))
    res = {}
    for vid, tr in spec:
        lo, hi = blob[pos:pos + n], blob[pos + n:pos + 2 * n]
        pos += 2 * n
        inter = bytearray(2 * n); inter[0::2] = lo; inter[1::2] = hi
        vals = list(struct.unpack('<%dh' % n, bytes(inter)))
        if tr == 1:
            vals = _undelta(vals)
        elif tr == 2:
            vals = [MISSING if x == MISSING else x + t for x, t in zip(vals, res['tmpf'])]
        res[VARS[vid]] = vals
    return fy, ny, res


# ── Roster ──────────────────────────────────────────────────────────────────
def load_roster(path_or_url):
    if re.match(r'https?://', path_or_url):
        with urllib.request.urlopen(path_or_url, timeout=60) as r:
            txt = r.read().decode('utf-8')
    else:
        txt = open(path_or_url, encoding='utf-8').read()
    rows = re.findall(r'\["([A-Z0-9]{3,4})","[^"]*","[A-Z]{2}",-?[\d.]+,-?[\d.]+,"([A-Z]{2})"\]', txt)
    exact = dict(re.findall(r'([A-Z0-9]{4}):\'([A-Za-z_/]+)\'', txt.split('STATION_TZ_EXACT', 1)[1].split('};', 1)[0])) \
        if 'STATION_TZ_EXACT' in txt else {}
    by_code = dict(re.findall(r'([A-Z]{2}):\'([A-Za-z_/]+)\'', txt.split('STATION_TZ_IANA', 1)[1].split('};', 1)[0])) \
        if 'STATION_TZ_IANA' in txt else {}
    stations, seen = [], set()
    for icao, code in rows:
        if icao in seen:
            continue
        seen.add(icao)
        stations.append({'id': icao, 'iana': exact.get(icao) or by_code.get(code) or 'America/New_York'})
    if len(stations) < 100:
        raise SystemExit(f'roster parse found only {len(stations)} stations; refusing to continue')
    return stations


# ── IEM fetching ────────────────────────────────────────────────────────────
class Transient(Exception):
    pass


class Throttled(Exception):
    def __init__(self, retry_after):
        super().__init__('HTTP 429'); self.retry_after = retry_after


THROTTLE_MAX = 12
STATS = {'requests': 0, 'throttled': 0, 'bytes': 0}


def iem_url(base, ids, sts, ets):
    q = '&'.join('station=' + s for s in ids) + '&' + '&'.join('data=' + f for f in FIELDS)
    return (f'{base}?{q}&sts={sts}&ets={ets}&tz=Etc%2FUTC&format=comma&latlon=no&elev=no'
            '&missing=M&trace=T&direct=no&report_type=3&report_type=4')


def http_get(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'contextclimate-hourly-archive (github.com/josht3113)',
                                               'Accept-Encoding': 'gzip'})
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            body = r.read()
            if r.headers.get('Content-Encoding') == 'gzip':
                body = gzip.decompress(body)
            return body.decode('utf-8', 'replace')
    except urllib.error.HTTPError as e:
        if e.code == 429:
            ra = e.headers.get('Retry-After') if e.headers else None
            raise Throttled(float(ra) if ra and ra.strip().isdigit() else 0)
        if e.code in (408, 425) or 500 <= e.code < 600:
            raise Transient(f'HTTP {e.code}')
        raise
    except (urllib.error.URLError, http.client.HTTPException, socket.timeout, TimeoutError, ConnectionError, OSError) as e:
        raise Transient(repr(e))


_last_request = [0.0]
def fetch_text(base, ids, sts, ets):
    """One request with retries; on repeated transient failure, split the station list.
    429 (rate limit) is waited out and never causes a split."""
    k = throttles = 0
    while True:
        wait = REQUEST_GAP_S - (time.time() - _last_request[0])
        if wait > 0:
            time.sleep(wait)
        _last_request[0] = time.time()
        STATS['requests'] += 1
        try:
            body = http_get(iem_url(base, ids, sts, ets))
            STATS['bytes'] += len(body)
            return body
        except Throttled as e:
            throttles += 1; STATS['throttled'] += 1
            if throttles > THROTTLE_MAX:
                raise SystemExit(f'IEM kept rate-limiting ({",".join(ids)} {sts}); stopping')
            time.sleep(max(e.retry_after, min(120, 15 * 2 ** min(throttles - 1, 3))))
            continue
        except Transient as e:
            if k == len(RETRY_WAITS):
                if len(ids) > 1:
                    h = len(ids) // 2
                    log(f'  splitting {len(ids)} stations after: {e}')
                    a = fetch_text(base, ids[:h], sts, ets)
                    b = fetch_text(base, ids[h:], sts, ets)
                    hdr_end = b.find('\n', b.find('station,')) + 1
                    return a + (b[hdr_end:] if b.find('station,') >= 0 else '')
                raise SystemExit(f'IEM request failed permanently ({",".join(ids)} {sts}..{ets}): {e}')
            log(f'  transient ({e}); retry in {RETRY_WAITS[k]} s')
            time.sleep(RETRY_WAITS[k])
            k += 1


def split_csv(text, wanted):
    """-> header index map, {station: [row parts]} (IEM drops the leading K for some ids)."""
    lines = text.split('\n')
    hdr = next((i for i, l in enumerate(lines) if l.startswith('station,')), -1)
    by = {s: [] for s in wanted}
    if hdr < 0:
        return None, by
    cols = lines[hdr].split(',')
    idx = {c: i for i, c in enumerate(cols)}
    alias = {}
    for s in wanted:
        alias[s] = s
        if len(s) == 4 and s[0] == 'K':
            alias[s[1:]] = s
    for l in lines[hdr + 1:]:
        if not l:
            continue
        p = l.split(',')
        sid = alias.get(p[0])
        if sid:
            by[sid].append(p)
    return idx, by


def build_decade(base, group, D, through):
    """Fetch UTC years D..min(D+9,through) (+ the first hours of D+10) for a station group."""
    y_hi = min(D + 9, through)
    y_lo = max(D, YEAR_FLOOR)
    if y_lo > y_hi:
        return {}
    per = {s['id']: [] for s in group}
    idx = None
    reqs = [(f'{Y}-01-01T00:00Z', f'{Y + 1}-01-01T00:00Z') for Y in range(y_lo, y_hi + 1)]
    reqs.append((f'{y_hi + 1}-01-01T00:00Z', f'{y_hi + 1}-01-01T14:00Z'))   # LST Dec 31 runs into UTC Jan 1
    for sts, ets in reqs:
        text = fetch_text(base, [s['id'] for s in group], sts, ets)
        i2, by = split_csv(text, per.keys())
        if i2 is not None:
            idx = i2
            for k, v in by.items():
                per[k].extend(v)
    out = {}
    for s in group:
        rows = per[s['id']]
        if not rows or idx is None:
            out[s['id']] = None
            continue
        rows.sort(key=lambda p: p[idx['valid']])
        years = parse_station(rows, idx, s['iana'], y_lo, y_hi)
        if all(v is None for v in years.values()):
            out[s['id']] = None
            continue
        out[s['id']] = (y_lo, [years[Y] for Y in range(y_lo, y_hi + 1)])
    return out


def write_decade(out_dir, sid, D, res):
    path = os.path.join(out_dir, 'hourly', sid, f'{D}.bin.gz')
    if res is None:
        if os.path.exists(path):
            os.remove(path)
        return False
    fy, years = res
    while years and years[0] is None:      # trim empty leading years
        years = years[1:]; fy += 1
    while years and years[-1] is None:
        years = years[:-1]
    if not years:
        return False
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'wb') as f:
        f.write(encode(fy, years))
    return True


def build_index(out_dir, through):
    stations = {}
    root = os.path.join(out_dir, 'hourly')
    for sid in sorted(os.listdir(root)) if os.path.isdir(root) else []:
        decs, first, last = [], None, None
        for fn in sorted(os.listdir(os.path.join(root, sid))):
            m = re.match(r'(\d{4})\.bin\.gz$', fn)
            if not m:
                continue
            with open(os.path.join(root, sid, fn), 'rb') as f:
                raw = gzip.decompress(f.read())
            fy, ny = struct.unpack('<HH', raw[6:10])
            decs.append(int(m.group(1)))
            first = fy if first is None else min(first, fy)
            last = fy + ny - 1 if last is None else max(last, fy + ny - 1)
        if decs:
            stations[sid] = {'first': first, 'last': last, 'decades': decs}
    idx = {'version': 2, 'built': dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%MZ'), 'through': through,
           'vars': VARS, 'missing': MISSING,
           'scale': {'tmpf': 100, 'dwpf': 100, 'feel': 100, 'vsby': 100, 'sky': 2, 'wind': 1, 'gust': 1, 'p01i': 100, 'wx': 1},
           'trace': TRACE_CODE, 'stations': stations}
    with open(os.path.join(out_dir, 'index.json'), 'w') as f:
        json.dump(idx, f, separators=(',', ':'))
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['full', 'incremental', 'stations', 'index'], required=True)
    ap.add_argument('--out', default='site')
    ap.add_argument('--stations', default='')
    ap.add_argument('--shard', default='0/1', help='k/n: build every n-th station starting at k')
    ap.add_argument('--through', type=int, default=0, help='last complete LST year (default: last calendar year)')
    ap.add_argument('--roster', default=ROSTER_URL)
    ap.add_argument('--iem', default=IEM_URL)
    ap.add_argument('--done-file', default='', help='resume log: "ICAO decade" lines already built')
    ap.add_argument('--max-minutes', type=float, default=0, help='stop cleanly (exit 3) after this long')
    a = ap.parse_args()

    through = a.through or dt.datetime.now(dt.timezone.utc).year - 1
    if a.mode == 'index':
        idx = build_index(a.out, through)
        log(f'index: {len(idx["stations"])} stations, through {through}')
        return

    roster = load_roster(a.roster)
    k, n = (int(x) for x in a.shard.split('/'))
    existing = {}
    ip = os.path.join(a.out, 'index.json')
    if os.path.exists(ip):
        existing = json.load(open(ip)).get('stations', {})

    plan = {}                                         # station id -> list of decades
    last_dec = through // 10 * 10
    all_decs = list(range(YEAR_FLOOR // 10 * 10, last_dec + 1, 10))
    wanted = set(x.strip().upper() for x in a.stations.split(',') if x.strip())
    for s in roster:
        if a.mode == 'stations':
            if s['id'] in wanted:
                plan[s['id']] = all_decs
        elif a.mode == 'full' or s['id'] not in existing:
            plan[s['id']] = all_decs
        else:
            plan[s['id']] = [last_dec]
    ids = [s['id'] for s in roster if s['id'] in plan][k::n]
    if a.mode == 'stations' and wanted - set(plan):
        log('not in roster:', ','.join(sorted(wanted - set(plan))))
    by_id = {s['id']: s for s in roster}
    log(f'mode {a.mode}, shard {a.shard}: {len(ids)} stations, through {through}')

    # Group stations that need the same decades, STATIONS_PER_REQUEST at a time.
    groups = {}
    for sid in ids:
        groups.setdefault(tuple(plan[sid]), []).append(by_id[sid])
    done = set()
    if a.done_file and os.path.exists(a.done_file):
        done = set(tuple(l.split()) for l in open(a.done_file) if l.strip())
        log(f'resuming: {len(done)} station-decades already built')
    written = 0
    t0 = time.time()
    for decs, members in groups.items():
        for g0 in range(0, len(members), STATIONS_PER_REQUEST):
            grp0 = members[g0:g0 + STATIONS_PER_REQUEST]
            for D in decs:
                grp = [s for s in grp0 if (s['id'], str(D)) not in done]
                if not grp:
                    continue
                if a.max_minutes and time.time() - t0 > a.max_minutes * 60:
                    log(f'time limit reached; stopping cleanly with {written} files written. Re-run to resume.')
                    sys.exit(3)
                td = time.time()
                res = build_decade(a.iem, grp, D, through)
                nw = sum(write_decade(a.out, s['id'], D, res.get(s['id'])) for s in grp)
                written += nw
                if nw:
                    log(f'  {D}s {grp[0]["id"]}..{grp[-1]["id"]}: {nw} files, {time.time() - td:.0f} s'
                        f' | total {STATS["requests"]} req, {STATS["bytes"] / 1e6:.0f} MB, {STATS["throttled"]} rate-limited')
                if a.done_file:
                    with open(a.done_file, 'a') as f:
                        f.writelines(f'{s["id"]} {D}\n' for s in grp)
            log(f'  {",".join(s["id"] for s in grp0)} done ({written} files, {time.time() - t0:.0f} s)')
    log(f'finished: {written} files in {time.time() - t0:.0f} s; {STATS}')
    if CLAMPED:
        log(f'values clamped to the int16 range (absurd source values): {CLAMPED}')


if __name__ == '__main__':
    main()
