#!/bin/bash
# --- Usage and defaults --------------------------------------------------------
print_usage() {
  printf "Usage: bash-example-input.sh --runpath /path/to/run --posctrl "barcode,ID barcode,ID" --negctrl "barcode,barcode"\n \
            options: \
               --help      Show this help message and exit.
               --runPath   Full path of the directory of the run folder. Required.
               --posCtrl   The positive controls for the run in the format of 'barcode,refID barcode,refID. Optional.'
               --negCtrl   The negative controls for the run in the format of 'barcode,barcode,barcode. Optional.'
  "
}

runPath=''
negCtrl=''
posCtrl=''

# --- Get arguments --------------------------------------------------------
OPTS=$(getopt -o '' -l "runPath:,posCtrl:,negCtrl:" --name "$0" -- "$@")

if [ $? -ne 0 ] ; then
        echo "Wrong input parameter!"; 1>&2
        print_usage
        exit 1;
fi

eval set -- "$OPTS"

# --- Process inputs --------------------------------------------------------
while true; do
    case "$1" in
        --runPath) runPath="$2"; shift 2 ;;
        --posCtrl) posCtrl="$2"; shift 2 ;;
        --negCtrl) negCtrl="$2"; shift 2 ;;
        --help) print_usage; exit 0 ;;
        --)        shift; break ;;
        *)         echo "Unexpected option: $1" >&2; exit 3 ;; 
    esac
done

# --- Ensure that --runPath is specified -------------------------------------
if [[ -z "$runPath" ]]; then
    echo "Error: --runPath is required" >&2
    print_usage
    exit 1
fi

# --- Pipeline goes here --------------------------------------------------------

echo "Run Path     | $runPath"
echo "Pos Controls | $posCtrl"
echo "Neg Controls | $negCtrl"