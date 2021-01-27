#!/bin/bash
mcs Tessel/*.cs CommonCode/MerStrings.cs -out:Tessel.exe
mono --aot -O=all,-shared Tessel.exe
