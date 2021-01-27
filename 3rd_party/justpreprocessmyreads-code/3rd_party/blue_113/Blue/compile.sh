#!/bin/bash
mcs Blue/*.cs CommonCode/MerStrings.cs -out:Blue.exe
mono --aot -O=all,-shared Blue.exe
