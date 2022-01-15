document.getElementById('re_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_re_cost').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_re_cost').setAttribute('hidden', true);
    }
});

document.getElementById('ligase_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_ligase_cost').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_ligase_cost').setAttribute('hidden', true);
    }
});