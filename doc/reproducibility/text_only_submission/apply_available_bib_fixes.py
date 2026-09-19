#!/usr/bin/env python3
"""Apply only unambiguous metadata additions from the supplied audit.

The input is never overwritten. Missing target entries stop the operation.
Unrecovered Shapely versions and original/translation choices are not guessed.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re

FIXES = {
    'kg_galeski2022lifshitz': ('pages', '7418'),
    'kg_shi2017lifshitz': ('pages', '14988'),
    'kg_li2019density': ('doi', '10.1134/S0081543819030076'),
}


def entries(text: str):
    cursor = 0
    pattern = re.compile(r'@([A-Za-z]+)\s*\{\s*([^,\s{}]+)\s*,')
    while match := pattern.search(text, cursor):
        opening = text.index('{', match.start())
        depth, quoted, escaped = 1, False, False
        index = opening + 1
        while index < len(text) and depth:
            char = text[index]
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == '"' and depth == 1:
                quoted = not quoted
            elif not quoted:
                depth += (char == '{') - (char == '}')
            index += 1
        if depth:
            raise ValueError('unterminated BibTeX entry')
        yield match.group(2), match.start(), match.end(), index - 1, index
        cursor = index


def field_ranges(body: str):
    depth, quoted, escaped, start = 0, False, False, 0
    for index, char in enumerate(body):
        if escaped:
            escaped = False
        elif char == '\\':
            escaped = True
        elif char == '"' and depth == 0:
            quoted = not quoted
        elif not quoted:
            depth += (char == '{') - (char == '}')
            if depth < 0:
                raise ValueError('invalid field nesting')
            if char == ',' and depth == 0:
                yield start, index
                start = index + 1
    if depth or quoted:
        raise ValueError('unterminated field value')
    if body[start:].strip():
        yield start, len(body)


def revise(text: str):
    found = {}
    for key, start, body_start, body_end, end in entries(text):
        if key in found:
            raise ValueError('duplicate bibliography key: ' + key)
        found[key] = (start, body_start, body_end, end)
    missing = sorted(set(FIXES) - found.keys())
    if missing:
        raise ValueError('target entries not found: ' + ', '.join(missing))
    replacements, changes = [], []
    for key, (field, value) in FIXES.items():
        _, left, right, _ = found[key]
        body = text[left:right]
        target = []
        for a, b in field_ranges(body):
            match = re.match(r'\s*([\w-]+)\s*=', body[a:b])
            if match and match.group(1).lower() == field:
                target.append((a, b))
        if len(target) > 1:
            raise ValueError('duplicate target field in ' + key)
        if target:
            a, b = target[0]
            old = body[a:b]
            prefix = re.match(r'\s*', old).group()
            suffix = re.search(r'\s*$', old).group()
            new = prefix + field + ' = {' + value + '}' + suffix
            replacements.append((left+a, left+b, new))
            changes.append({'key': key, 'field': field, 'old': old.strip(), 'new': value})
        else:
            prefix = '' if not body.strip() or body.rstrip().endswith(',') else ','
            replacements.append((right, right, prefix + '\n  ' + field + ' = {' + value + '}\n'))
            changes.append({'key': key, 'field': field, 'old': None, 'new': value})
    for left, right, value in sorted(replacements, reverse=True):
        text = text[:left] + value + text[right:]
    return text, changes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.input.resolve() == args.output.resolve():
        parser.error('use a new output path; the input will not be overwritten')
    original = args.input.read_text(encoding='utf-8')
    updated, changes = revise(original)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(updated, encoding='utf-8')
    report = {'input': str(args.input), 'output': str(args.output),
              'input_sha256': sha256(original.encode()).hexdigest(),
              'output_sha256': sha256(updated.encode()).hexdigest(),
              'changes': changes, 'unresolved_choices_modified': False}
    args.output.with_suffix(args.output.suffix + '.changes.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
