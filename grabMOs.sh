#!/bin/bash


function getMos {
    cat $1 | sed -n '/Index Energy  Occupation Coefficients/,/--/p' > orbitals.txt
}