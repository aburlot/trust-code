#!/bin/bash
echo "Usage: gdb.sh [-valgrind] \${exec}"

valgrind=0
if [[ "${1}" = -valgrind ]]; then
    valgrind=1
    shift
fi

directory=""
if [[ "${1}" != "" ]] && [[ -f "${1}" ]]; then
    #Le lancement prend desormais en compte un repertoire d'include pour gdb 6.x
    exec=`basename $1`
    MonoDir=$TRUST_ROOT/MonoDir${exec#TRUST}
    if [[ -d "${MonoDir}/src" ]]; then
        directory="--directory=${MonoDir}/src --directory=${TRUST_ROOT}/include"
    fi

    # Pour prendre en compte l'atelier on ajoute aussi le repertoire du binaire
    directory="${directory} --directory=$(dirname "${1}")"
fi

if [[ "${valgrind}" == 1 ]]; then
    valgrind --vgdb=yes --vgdb-error=0 --leak-check=full --track-origins=yes --show-reachable=yes --suppressions="${TRUST_ROOT}/Outils/valgrind/suppressions" --num-callers=15 "${@}" &
    pid=$!
    sleep 1
    "${Xterm}" -e /usr/bin/gdb "${directory}" "${@}" -ex "target remote | $TRUST_ROOT/exec/valgrind/lib/valgrind/../../bin/vgdb --pid=${pid}" -ex "echo Please wait a few seconds... The valgrind process is slow to start.\n " -ex continue
else
    # Surcharge de gdb pour pouvoir afficher des tableaux
    # Lancement de valgrind
    options=${TRUST_ROOT}/bin/KSH/gdb.options

    echo "---------------------------------------------------------------------"
    grep "To print" <"${options}"
    echo "# To debug dynamically-linked executable, use: set auto-solib-add off"
    echo "Then load the library you want to debug: share name_of_library.so"
    echo "---------------------------------------------------------------------"
    echo "Lancement de: /usr/bin/gdb ${directory} -x ${options} ${@}"

    /usr/bin/gdb ${directory} -x ${options} ${@}
fi
