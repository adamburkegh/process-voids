'''
Quick sanity check for a result CSV - the "did this run go the way I
expected" question, answered once instead of as a fresh python -c
snippet each time: row count, status breakdown, and which logs/combos/
degradation dims the CSV actually covers. A column that isn't present
(eg a *_timings.csv has no 'combo'/'degradation_dim') is skipped, not
an error - not every result CSV shares the same schema.
'''

import argparse

import pandas as pd

_LIST_COLUMNS = ['log', 'combo', 'degradation_dim']


def summarize(csv_path):
    df = pd.read_csv(csv_path)
    summary = {'rows': len(df)}
    summary['status'] = df['status'].value_counts().to_dict() if 'status' in df.columns else None
    for col in _LIST_COLUMNS:
        key = 'degradation_dims' if col == 'degradation_dim' else f'{col}s'
        summary[key] = sorted(df[col].unique()) if col in df.columns else None
    return summary


def format_summary(summary):
    lines = [f"rows: {summary['rows']}"]
    for key, value in summary.items():
        if key != 'rows' and value is not None:
            lines.append(f'{key}: {value}')
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='Quick sanity summary of a result CSV - row count, status, coverage.')
    parser.add_argument('csv', help='Result CSV path')
    args = parser.parse_args()
    print(format_summary(summarize(args.csv)))


if __name__ == '__main__':
    main()
