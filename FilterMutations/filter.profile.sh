PYDIR=/home/jocostello/shared/LG3_Pipeline
if ! echo ${PYTHONPATH} | egrep "(^|:)${PYDIR}($|:)" >/dev/null ; then
	export PYTHONPATH=${PYTHONPATH}:$PYDIR
fi

DIR=${PYDIR}/FilterMutations
if [ -d "${DIR}" ]; then
    if ! echo ${PATH} | egrep "(^|:)$DIR($|:)" >/dev/null ; then
        export PATH=$PATH:$DIR
    fi
fi

