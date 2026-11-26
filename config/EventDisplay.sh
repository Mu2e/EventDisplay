#!/bin/bash

# =======================
# Mu2e Event Display Script
# =======================


usage() {
  echo "Usage: $0
  [ --DATASET ]
  [ --RUN ]
  [ --SUBRUN  ]
  [ --EVENT ]
  [ --GEOM ]
  ./EventDisplay.sh --run 1201 --subrun 34 --event 15028 --geom 2025 --dataset mcs.mu2e.ensembleMDS2cMix1BBTriggered.MDC2020ba_best_v1_3.art
  or ./EventDisplay.sh  --run 1201 --subrun 476 --geom 2025 --event 1  --dataset nts.mu2e.ensembleMDS2cMix1BBTriggered.MDC2020ba_best_v1_3_v06_06_00.001201_00000476.root
  " 1>&2
}



# Defaults, overridden by the command-line options below
RUN=0
SUBRUN=0
EVENT=0
DATASET=""
GEOM=2025

# Function: Exit with error.
exit_abnormal() {
  usage
  exit 1
}


# Loop: Get the next option;
while getopts ":-:" options; do
  case "${options}" in
    -)
      case "${OPTARG}" in
        dataset)
          DATASET=${!OPTIND} OPTIND=$(( $OPTIND + 1 ))
          ;;
        event)
          EVENT=${!OPTIND} OPTIND=$(( $OPTIND + 1 ))
          ;;
        subrun)
          SUBRUN=${!OPTIND} OPTIND=$(( $OPTIND + 1 ))
          ;;
        run)
          RUN=${!OPTIND} OPTIND=$(( $OPTIND + 1 ))
          ;;
        geom)
          GEOM=${!OPTIND} OPTIND=$(( $OPTIND + 1 ))
          ;;
      esac;;
    :)                                    # If expected argument omitted:
      echo "Error: -${OPTARG} requires an argument."
      exit_abnormal                       # Exit abnormally.
      ;;
    *)                                    # If unknown (any other) option:
      exit_abnormal                       # Exit abnormally.
      ;;
  esac
done

echo "========== Launching Mu2e/EventDisplay ==============="

# Command 1: pickEvent
echo "Extracting the data-set..."

# Construct the EVENT identifier string in the format "RUN/SUBRUN/EVENT"
RUN_SUBRUN_EVENT="${RUN}/${SUBRUN}/${EVENT}"

first_three_chars="${DATASET:0:3}"

echo "filetype is ${first_three_chars}"

if [ "$first_three_chars" == "nts" ]; then
  FILE=$(samweb file-lineage parents ${DATASET})
  echo "file is ${FILE}"
  NEWNAME=${FILE%.*}
  NEWNAME=${NEWNAME%.*}
  DATASET=${NEWNAME}.art
  echo "new name is ${DATASET}"
fi

# Execute the pickEvent command
# If the command fails, the script will terminate with an error message
pickEvent -e -v "${DATASET}" "${RUN_SUBRUN_EVENT}"
if [ $? -ne 0 ]; then
    echo "Error: pickEvent command failed."
    exit 1
fi

# Command 2: mu2e
echo "Data-set has been extracted, now RUNning display."


# Construct the art file name
art_file_name="${DATASET}_${RUN}_${SUBRUN}_${EVENT}.art"

# Execute the mu2e command
# If the command fails, the script will terminate with an error message
mu2e -c EventDisplay/examples/nominal_MDC${GEOM}.fcl ${art_file_name}
if [ $? -ne 0 ]; then
    echo "Error: mu2e command failed."
    exit 1
fi

echo "Display finished."

