"""Minimal xlsx sheet reader (cell values only) that ignores drawings.

openpyxl fails on the SiVIRA workbook's broken drawing reference, so build_national.py reads it with
this instead. read(path, sheet) returns one dict per row, keyed by column letter.
Usage on its own: python3 -I xlsx_read.py file.xlsx sheetname [n_rows]
"""
import zipfile, re, sys, xml.etree.ElementTree as ET

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}


def read(path, sheet):
    z = zipfile.ZipFile(path)
    wb = ET.fromstring(z.read('xl/workbook.xml'))
    rels = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    rid = {r.get('Id'): r.get('Target') for r in rels}
    target = None
    for s in wb.find('m:sheets', NS):
        if s.get('name') == sheet:
            target = rid[s.get('{%s}id' % NS['r'])]
    target = target.lstrip('/')
    target = target if target.startswith('xl/') else 'xl/' + target
    ss = []
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', NS):
            ss.append(''.join(t.text or '' for t in si.iter('{%s}t' % NS['m'])))
    out = []
    for row in ET.fromstring(z.read(target)).iter('{%s}row' % NS['m']):
        vals = {}
        for c in row.findall('m:c', NS):
            col = re.match(r'[A-Z]+', c.get('r')).group(0)
            v = c.find('m:v', NS)
            if v is None:
                t = c.find('m:is/m:t', NS)
                vals[col] = t.text if t is not None else None
                continue
            vals[col] = ss[int(v.text)] if c.get('t') == 's' else v.text
        out.append(vals)
    return out


if __name__ == '__main__':
    for r in read(sys.argv[1], sys.argv[2])[:int(sys.argv[3]) if len(sys.argv) > 3 else 5]:
        print(r)
