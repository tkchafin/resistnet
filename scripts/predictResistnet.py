#!/usr/bin/env python3
"""Apply a fitted ResistNet model to the reaches of another network."""
import argparse
import sys

import geopandas as gpd

import resistnet.predict as pr


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Apply a fitted model ($out.modelSpec.json from runResistnet.py) '
            'to the reaches of another network, e.g. another watershed. '
            'Variables are rescaled and transformed with the ranges of the '
            'network the model was fitted to; values outside those ranges '
            'are clamped to them.'
        )
    )
    parser.add_argument(
        '-s', '--shapefile', required=True,
        help='Reaches (shapefile or GeoPackage) with the model\'s variables'
    )
    parser.add_argument(
        '-m', '--model_spec', required=True,
        help='Model specification ($out.modelSpec.json)'
    )
    parser.add_argument(
        '-n', '--model', type=int, default=1,
        help=('Model to apply, by its rank in the Hall of Fame '
              '(default: 1, the best)')
    )
    parser.add_argument(
        '--reachid_col', default=None,
        help='Reach ID attribute (default: as when the model was fitted)'
    )
    parser.add_argument(
        '--length_col', default=None,
        help='Reach length attribute (default: as when the model was fitted)'
    )
    parser.add_argument(
        '--next_down_col', default='NEXT_DOWN',
        help=('SAMC: attribute giving the next downstream reach ID '
              '(default: NEXT_DOWN)')
    )
    parser.add_argument(
        '-o', '--out', default='predicted',
        help='Output prefix (default: predicted)'
    )
    args = parser.parse_args()

    spec = pr.load_spec(args.model_spec)
    if not 1 <= args.model <= len(spec["models"]):
        sys.exit(f"--model must be from 1 to {len(spec['models'])}")
    reaches = gpd.read_file(args.shapefile)
    table, moves, clamped = pr.predict_reaches(
        reaches, spec, args.model - 1, reachid_col=args.reachid_col,
        next_down_col=args.next_down_col, length_col=args.length_col)

    kind = spec["models"][args.model - 1]["type"]
    for var, share in clamped.items():
        if share > 0 and var != "composite":
            print(f"WARNING: {share:.1%} of reaches have {var} outside the "
                  "range of the network the model was fitted to; their "
                  "values were clamped to that range")
    if clamped.get("composite", 0) > 0:
        unit = "reaches" if kind == "standard" else "moves between reaches"
        print(f"WARNING: on {clamped['composite']:.1%} of {unit}, the "
              "variables combine to a value outside the range of the network "
              "the model was fitted to; it was clamped to that range, so "
              "those values don't differ from each other")

    for col in table.columns[1:]:
        reaches[col] = table[col].to_numpy()
    reaches.to_file(f"{args.out}.predicted.gpkg", driver="GPKG")
    table.to_csv(f"{args.out}.predicted.tsv", sep="\t", index=False,
                 na_rep="NA")
    if moves is not None:
        moves.to_csv(f"{args.out}.predictedMoves.tsv", sep="\t", index=False)
    print(f"Applied model {args.model} ({kind}) to {len(reaches)} reaches; "
          f"outputs: {args.out}.predicted.*")


if __name__ == "__main__":
    main()
