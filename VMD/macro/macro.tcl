# POPC head and tail, the head includes also O22 O21 to allow MSMS to easily locate the surface
atomselect macro popchead {
    (lipid and not name "C[23].+" "H[0-9]+[RSXYTZ]" H91 H101)
}
atomselect macro popctail {
    (lipid and name "C[23].+" "H[0-9]+[RSXYTZ]" H91 H101)
}
# KcsA filter (TVGYG sequence)
atomselect macro kcsafilter {
    (protein and (( resid 74 to 79 and backbone) or (resid 75 and noh)))
}


