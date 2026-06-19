#! /bin/bash

# Runs a repo-local ltsep wrapper so batch jobs do not depend on upstream
# getPathDict.py calling os.getlogin() on worker nodes.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR}"
LTSEP_FIELDS_SCRIPT="${REPO_ROOT}/farm_env/print_ltsep_path_fields.py"
if [[ ! -f "${LTSEP_FIELDS_SCRIPT}" ]]; then
    REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
    LTSEP_FIELDS_SCRIPT="${REPO_ROOT}/farm_env/print_ltsep_path_fields.py"
fi
PATHFILE_INFO="$(python3 "${LTSEP_FIELDS_SCRIPT}" "${REPO_ROOT}")"
path_rc=$?
if [[ "${path_rc}" -ne 0 || -z "${PATHFILE_INFO}" ]]; then
    echo "ERROR: failed to resolve ltsep paths via ${LTSEP_FIELDS_SCRIPT}" >&2
    exit 1
fi

# Split the string we get to individual variables, easier for printing and use later
VOLATILEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
ANALYSISPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f2`
HCANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f3`
REPLAYPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f4`
UTILPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f5`
PACKAGEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f6`
OUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f7`
ROOTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f8`
SKIMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f9`
REPORTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f10`
CUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f11`
PARAMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f12`
SCRIPTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f13`
ANATYPE=`echo ${PATHFILE_INFO} | cut -d ','  -f14`
USER=`echo ${PATHFILE_INFO} | cut -d ','  -f15`
HOST=`echo ${PATHFILE_INFO} | cut -d ','  -f16`
SIMCPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f17`
LTANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f18`

sync_symlink() {
    local link_path="$1"
    local target_path="$2"
    local label="$3"
    local current_target=""

    if [ -L "${link_path}" ]; then
        current_target=$(readlink "${link_path}")
        if [ "${current_target}" = "${target_path}" ] && [ -e "${link_path}" ]; then
            echo "${label} sym link already exists and not broken"
            echo "             ${link_path}-->${target_path}"
            echo
            echo
            return 0
        fi

        echo "${label} sym link exists but needs updating, replacing"
        rm "${link_path}"
    elif [ -e "${link_path}" ]; then
        echo "ERROR: ${link_path} exists and is not a symlink"
        return 1
    fi

    ln -s "${target_path}" "${link_path}"
}


# Flag definitions (flags: h, a, o, s)
while getopts 'h' flag; do
    case "${flag}" in
        h) 
        echo "--------------------------------------------------------------"
        echo "./set_SymLinks.sh {variable arguments, see help}"
	echo
        echo "Description: Setup symlinks across lt_analysis and simc_gfortran"
        echo "--------------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."
	echo "    ParticleType=arg1"
        echo "    -h, help"
        exit 0
        ;;
        *) print_usage
        exit 1 ;;
    esac
done

ParticleType=$1
BACKGROUND_SAMPLE_DIR="${LTANAPATH}/background_samples"
BACKGROUND_CONFIG_FILE="${BACKGROUND_SAMPLE_DIR}/background_samples.conf"
BACKGROUND_OUTPUT_DIR="${VOLATILEPATH}/OUTPUT/Analysis/${ANATYPE}LT/background_samples/"
BACKGROUND_OUTPUT_LINK="${BACKGROUND_SAMPLE_DIR}/OUTPUTS"
BACKGROUND_SIMC_LINK="${BACKGROUND_SAMPLE_DIR}/simc_gfortran"
BACKGROUND_SIMC_PATH="${BACKGROUND_SIMC_PATH:-}"

if [ -f "${BACKGROUND_CONFIG_FILE}" ]; then
    # shellcheck disable=SC1090
    source "${BACKGROUND_CONFIG_FILE}"
fi

mkdir -p "${BACKGROUND_OUTPUT_DIR}"

# Check symlinks and create/fix if bad
if [ ! -L "${LTANAPATH}/simc_gfortran" ]; then
    ln -s "${SIMCPATH}" "${LTANAPATH}/simc_gfortran"
elif [ -L "${LTANAPATH}/simc_gfortran" ]; then
    if [ ! -e "${LTANAPATH}/simc_gfortran" ]; then
	echo "${LTANAPATH}/simc_gfortran sym link exits but is broken, replacing"
	rm "${LTANAPATH}/simc_gfortran"
	ln -s "${SIMCPATH}" "${LTANAPATH}/simc_gfortran"
    else 
	echo "${LTANAPATH}/simc_gfortran sym link already exists and not broken"
	echo "             ${LTANAPATH}/simc_gfortran-->${SIMCPATH}"
	echo
	echo
    fi
fi

if [ ! -L "${LTANAPATH}/OUTPUT" ]; then
    ln -s "${VOLATILEPATH}/OUTPUT/" "${LTANAPATH}/OUTPUT"
elif [ -L "${LTANAPATH}/OUTPUT" ]; then
    if [ ! -e "${LTANAPATH}/OUTPUT" ]; then
	echo "${LTANAPATH}/OUTPUT sym link exits but is broken, replacing"
	rm "${LTANAPATH}/OUTPUT"
	ln -s "${VOLATILEPATH}/OUTPUT/" "${LTANAPATH}/OUTPUT"
    else 
	echo "${LTANAPATH}/OUTPUT sym link already exists and not broken"
	echo "             ${LTANAPATH}/OUTPUT-->${VOLATILEPATH}/OUTPUT/"
	echo
	echo
    fi
fi

if [ ! -L "${LTANAPATH}/log" ]; then
    ln -s "${VOLATILEPATH}/log/" "${LTANAPATH}/log"
elif [ -L "${LTANAPATH}/log" ]; then
    if [ ! -e "${LTANAPATH}/log" ]; then
	echo "${LTANAPATH}/log sym link exits but is broken, replacing"
	rm "${LTANAPATH}/log"
	ln -s "${VOLATILEPATH}/log/" "${LTANAPATH}/log"
    else 
	echo "${LTANAPATH}/log sym link already exists and not broken"
	echo "             ${LTANAPATH}/log-->${VOLATILEPATH}/log/"
	echo
	echo
    fi
fi

if [ ! -L "${LTANAPATH}/src/data" ]; then
    ln -s "${VOLATILEPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "${LTANAPATH}/src/data"
elif [ -L "${LTANAPATH}/src/data" ]; then
    if [ ! -e "${LTANAPATH}/src/data" ]; then
	echo "${LTANAPATH}/src/data sym link exits but is broken, replacing"
	rm "${LTANAPATH}/src/data"
	ln -s "${VOLATILEPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "${LTANAPATH}/src/data"
    else 
	echo "${LTANAPATH}/src/data sym link already exists and not broken"
	echo "             ${LTANAPATH}/src/data-->${VOLATILEPATH}/OUTPUT/Analysis/${ANATYPE}LT/"
	echo
	echo
    fi
fi

if [ ! -L "${LTANAPATH}/simc_gfortran/OUTPUTS" ]; then
    ln -s "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "${LTANAPATH}/simc_gfortran/OUTPUTS"
elif [ -L "${LTANAPATH}/simc_gfortran/OUTPUTS" ]; then
    if [ ! -e "${LTANAPATH}/simc_gfortran/OUTPUTS" ]; then
	echo "${LTANAPATH}/simc_gfortran/OUTPUTS sym link exits but is broken, replacing"
	rm "${LTANAPATH}/simc_gfortran/OUTPUTS"
	ln -s "${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/" "${LTANAPATH}/simc_gfortran/OUTPUTS"
    else 
	echo "${LTANAPATH}/simc_gfortran/OUTPUTS sym link already exists and not broken"
	echo "             ${LTANAPATH}/simc_gfortran/OUTPUTS-->${LTANAPATH}/OUTPUT/Analysis/${ANATYPE}LT/"
	echo
	echo
    fi
fi

if [ ! -L "${LTANAPATH}/simc_gfortran/input" ]; then
    ln -s "${LTANAPATH}/input/${ParticleType}/" "${LTANAPATH}/simc_gfortran/input"
elif [ -L "${LTANAPATH}/simc_gfortran/input" ]; then
    if [ ! -e "${LTANAPATH}/simc_gfortran/input" ]; then
	echo "${LTANAPATH}/simc_gfortran/input sym link exits but is broken, replacing"
	rm "${LTANAPATH}/simc_gfortran/input"
	ln -s "${LTANAPATH}/input/${ParticleType}/" "${LTANAPATH}/simc_gfortran/input"
    else 
	echo "${LTANAPATH}/simc_gfortran/input sym link already exists and not broken"
	echo "             ${LTANAPATH}/simc_gfortran/input-->${LTANAPATH}/input/${ParticleType}/"
	echo
	echo
    fi
fi

if [ ! -L "${LTANAPATH}/simc_gfortran/worksim" ]; then
    ln -s "${VOLATILEPATH}/worksim/" "${LTANAPATH}/simc_gfortran/worksim"
elif [ -L "${LTANAPATH}/simc_gfortran/worksim" ]; then
    if [ ! -e "${LTANAPATH}/simc_gfortran/worksim" ]; then
	echo "${LTANAPATH}/simc_gfortran/worksim sym link exits but is broken, replacing"
	rm "${LTANAPATH}/simc_gfortran/worksim"
	ln -s "${VOLATILEPATH}/worksim/" "${LTANAPATH}/simc_gfortran/worksim"
    else 
	echo "${LTANAPATH}/simc_gfortran/worksim sym link already exists and not broken"
	echo "             ${LTANAPATH}/simc_gfortran/worksim-->${VOLATILEPATH}/worksim/"
	echo
	echo
	    fi
fi

if [ -n "${BACKGROUND_SIMC_PATH}" ]; then
    sync_symlink "${BACKGROUND_SIMC_LINK}" "${BACKGROUND_SIMC_PATH}" "${BACKGROUND_SIMC_LINK}" || exit 1
fi

sync_symlink "${BACKGROUND_OUTPUT_LINK}" "${BACKGROUND_OUTPUT_DIR}" "${BACKGROUND_OUTPUT_LINK}" || exit 1

if [ -L "${BACKGROUND_SIMC_LINK}" ] && [ -e "${BACKGROUND_SIMC_LINK}" ]; then
    sync_symlink "${BACKGROUND_SIMC_LINK}/worksim" "${VOLATILEPATH}/worksim/" "${BACKGROUND_SIMC_LINK}/worksim" || exit 1
else
    echo "Background SIMC symlink is not configured yet"
    echo "             Set BACKGROUND_SIMC_PATH in ${BACKGROUND_CONFIG_FILE} to enable it"
    echo
    echo
fi

echo 
echo "Directories and sym links for ${LTANAPATH} now setup"

exit 0
