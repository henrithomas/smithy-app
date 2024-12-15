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

document.getElementById('ecori_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_EcoRI_cost').removeAttribute('hidden');
        document.getElementById('id_EcoRI_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_EcoRI_cost').setAttribute('hidden', true);
        document.getElementById('id_EcoRI_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('xbai_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_XbaI_cost').removeAttribute('hidden');
        document.getElementById('id_XbaI_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_XbaI_cost').setAttribute('hidden', true);
        document.getElementById('id_XbaI_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('spei_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_SpeI_cost').removeAttribute('hidden');
        document.getElementById('id_SpeI_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_SpeI_cost').setAttribute('hidden', true);
        document.getElementById('id_SpeI_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('psti_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_PstI_cost').removeAttribute('hidden');
        document.getElementById('id_PstI_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_PstI_cost').setAttribute('hidden', true);
        document.getElementById('id_PstI_n_reacts').setAttribute('hidden', true);
    }
});