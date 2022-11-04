document.getElementById('mmix_check').addEventListener('click', function()
{
    if(this.checked)
    {
        document.getElementById('mmix_col').removeAttribute('hidden');
        document.getElementById('pcr_poly_col').setAttribute('hidden', true);

    }
    else
    {
        document.getElementById('mmix_col').setAttribute('hidden', true);
        document.getElementById('pcr_poly_col').removeAttribute('hidden');
    }
});

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