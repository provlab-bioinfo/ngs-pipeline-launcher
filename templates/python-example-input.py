import argparse, warnings

def posCtrlType(s):
    try:
        si = s.split(',')
        si[0] = int(si[0])
        return tuple(si)
    except:
        raise argparse.ArgumentTypeError("Positive controls must be given divided by commas and space, dot, or semicolon e.g.: '94,21-117-036013 95,21-129-012673'")

parser = argparse.ArgumentParser()

parser.add_argument("--runPath", required=True, type=str, help="Full path of the directory of the run folder.")
parser.add_argument("--posCtrl", type=posCtrlType, nargs='+', help="The positive controls for the run in the format of 'barcode,refID barcode,refID. Optional.'")
parser.add_argument("--negCtrl", type=str, help="The negative controls for the run in the format of 'barcode,barcode,barcode'. Optional.")

args, unknown = parser.parse_known_args()

if unknown:
    bad_args = ', '.join(unknown)    
    warnings.warn(f"Ignoring unrecognized arguments: {bad_args}", UserWarning)

print(f"Run Path     | {args.runPath}")
print(f"Pos Controls | {args.posCtrl}")
print(f"Neg Controls | {args.negCtrl}")