#!/bin/bash


function getScfOrbitals {
    cat "$1" | sed -n '/Index Energy  Occupation Coefficients/,/--/p' > orbitals.txt
}

function getRasOrbitals {
    sed -n '/INDEX  ENERGY  OCCUPATION COEFFICIENTS .../,/^--$/p' "$1" > orbitals.txt
}