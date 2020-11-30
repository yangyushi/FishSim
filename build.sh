export prefix="/usr/bin/local"
cd src
make clean
make
if [[ $1 == "test" ]]; then
    cd ../script
    python3 test_csimulate.py
fi
cd ..
