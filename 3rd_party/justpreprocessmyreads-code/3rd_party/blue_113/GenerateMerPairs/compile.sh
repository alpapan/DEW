#!/bin/bash
mcs GenerateMerPairs/*.cs Tessel/MerCollections.cs CommonCode/MerStrings.cs -out:GenerateMerPairs.exe
mono --aot -O=all,-shared GenerateMerPairs.exe
