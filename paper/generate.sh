docker run --rm \
    --volume $PWD/paper:/data \
    --env JOURNAL=joss \
    openjournals/inara

# --user $(id -u):$(id -g) \
