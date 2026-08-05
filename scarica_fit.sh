#!/bin/sh
# Scarica i fit delle repliche dalla release GitHub e li scompatta in output/.
# Richiede la GitHub CLI autenticata (https://cli.github.com, poi `gh auth login`).
# In alternativa: scaricare a mano i .tar dalla pagina della release
# e scompattarli con `tar xf <file>.tar -C output`.
set -e
TAG="${1:-fit-v1}"
mkdir -p output
gh release download "$TAG" --pattern '*.tar' --dir output
for f in output/*.tar; do
    tar xf "$f" -C output && rm "$f"
done
echo "Fatto: fit in output/fit/, risultati derivati in output/."
