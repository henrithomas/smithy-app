document.getElementById('pcr_poly_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_pcr_polymerase_cost').removeAttribute('hidden');
        document.getElementById('id_pcr_polymerase_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_pcr_polymerase_cost').setAttribute('hidden', true);
        document.getElementById('id_pcr_polymerase_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('exo_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_exonuclease_cost').removeAttribute('hidden');
        document.getElementById('id_exonuclease_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_exonuclease_cost').setAttribute('hidden', true);
        document.getElementById('id_exonuclease_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('ligase_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_ligase_cost').removeAttribute('hidden');
        document.getElementById('id_ligase_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_ligase_cost').setAttribute('hidden', true);
        document.getElementById('id_ligase_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('poly_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_polymerase_cost').removeAttribute('hidden');
        document.getElementById('id_polymerase_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_polymerase_cost').setAttribute('hidden', true);
        document.getElementById('id_polymerase_n_reacts').setAttribute('hidden', true);
    }
});