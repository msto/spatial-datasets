#!/usr/bin/env awk

BEGIN {FS=OFS="\t"}

{
    if (NR==1) {
        $9="Cluster";
    } else if ($6 <= 0) {
        $9="Unassigned";
    } else if ($6 <= 5) {
        $9="Malignant";
    } else if ($6 <= 9) {
        $9="Fibroblast";
    } else if ($6 <= 13) {
        $9="Macrophage";
    } else if ($6 <= 15) {
        $9="DC";
    } else if ($6 == 16) {
        $9="B-cells";
    } else if ($6 == 17) {
        $9="T-cells";
    } else if ($6 == 18) {
        $9="Erythrocytes";
    } else {
        $9="Unassigned"
    }

    print
}
