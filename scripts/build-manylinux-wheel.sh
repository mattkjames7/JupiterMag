#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: scripts/build-manylinux-wheel.sh [python-version]

Build and repair a wheel in the same manylinux container used by CI.

Examples:
  scripts/build-manylinux-wheel.sh
  scripts/build-manylinux-wheel.sh 3.14
  scripts/build-manylinux-wheel.sh 3.12

Outputs:
  dist/*.whl
  wheelhouse/*.whl
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

PYVER="${1:-3.14}"

case "${PYVER}" in
  3.10) PYTAG="cp310-cp310" ;;
  3.11) PYTAG="cp311-cp311" ;;
  3.12) PYTAG="cp312-cp312" ;;
  3.13) PYTAG="cp313-cp313" ;;
  3.14) PYTAG="cp314-cp314" ;;
  *)
    echo "Unsupported python-version: ${PYVER}" >&2
    usage
    exit 1
    ;;
esac

if command -v git >/dev/null 2>&1; then
  REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
else
  REPO_ROOT="$(pwd)"
fi

IMAGE="quay.io/pypa/manylinux_2_34:latest"

echo "Building manylinux wheel for Python ${PYVER} (${PYTAG})"
echo "Repository root: ${REPO_ROOT}"

rm -rf "${REPO_ROOT}/dist" "${REPO_ROOT}/wheelhouse"
mkdir -p "${REPO_ROOT}/dist" "${REPO_ROOT}/wheelhouse"

docker run --rm \
  -e PYTAG="${PYTAG}" \
  -e HOST_UID="$(id -u)" \
  -e HOST_GID="$(id -g)" \
  -v "${REPO_ROOT}:/io" \
  -w /io \
  "${IMAGE}" \
  /bin/bash -lc '
    set -euo pipefail
    if ! command -v cmake >/dev/null 2>&1; then
      dnf install -y cmake gcc gcc-c++ make ninja-build
    fi
    "/opt/python/${PYTAG}/bin/python" -m pip install --upgrade pip setuptools wheel build auditwheel
    "/opt/python/${PYTAG}/bin/python" -m build --wheel -o dist
    auditwheel repair dist/*.whl -w wheelhouse
    chown -R "${HOST_UID}:${HOST_GID}" /io/dist /io/wheelhouse || true
  '

echo "Done. Built files:"
ls -1 "${REPO_ROOT}/wheelhouse"/*.whl