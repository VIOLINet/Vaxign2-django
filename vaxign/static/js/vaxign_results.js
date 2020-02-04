function chk_opt(elementID1, elementID2) {
    var elementSel1 = document.getElementById(elementID1);
    var elementSel2 = document.getElementById(elementID2);
    for(j=0; j<elementSel1.options.length; j++) {
        if (elementSel1.options[j].selected) {
            elementSel2.options[j].selected = false;
        }
    }
}

function show_selected_no(elementID1, elementID2) {
    var elementSel1 = document.getElementById(elementID1);
    var elementSel2 = document.getElementById(elementID2);
    var count=0;
    for(j=0; j<elementSel1.options.length; j++) {
        if (elementSel1.options[j].selected) {
            count++;
        }
    }
    elementSel2.innerHTML=count;
}