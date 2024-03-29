document.getElementById('id_gibson').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('gibson_input').removeAttribute('hidden');
    }
    else {
        document.getElementById('gibson_input').setAttribute('hidden', true);
    }
});

document.getElementById('id_goldengate').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('goldengate_input').removeAttribute('hidden');
    }
    else {
        document.getElementById('goldengate_input').setAttribute('hidden', true);
    }
});

document.getElementById('id_pcr').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('pcr_input').removeAttribute('hidden');
    }
    else {
        document.getElementById('pcr_input').setAttribute('hidden', true);
    }
});

document.getElementById('id_slic').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('slic_input').removeAttribute('hidden');
    }
    else {
        document.getElementById('slic_input').setAttribute('hidden', true);
    }
});

document.getElementById('id_biobricks').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('biobricks_input').removeAttribute('hidden');
    }
    else {
        document.getElementById('biobricks_input').setAttribute('hidden', true);
    }
});

document.getElementById('gib_mmix_check').addEventListener('click', function()
{
    if(this.checked)
    {
        document.getElementById('gib_mmix_col').removeAttribute('hidden');
        document.getElementById('gib_exo_col').setAttribute('hidden', true);
        document.getElementById('gib_lig_col').setAttribute('hidden', true);
        document.getElementById('gib_poly_col').setAttribute('hidden', true);

    }
    else
    {
        document.getElementById('gib_mmix_col').setAttribute('hidden', true);
        document.getElementById('gib_exo_col').removeAttribute('hidden');
        document.getElementById('gib_lig_col').removeAttribute('hidden');
        document.getElementById('gib_poly_col').removeAttribute('hidden');
    }
});

document.getElementById('gg_mmix_check').addEventListener('click', function()
{
    if(this.checked)
    {
        document.getElementById('gg_mmix_col').removeAttribute('hidden');
        document.getElementById('gg_re_col').setAttribute('hidden', true);
        document.getElementById('gg_lig_col').setAttribute('hidden', true);

    }
    else
    {
        document.getElementById('gg_mmix_col').setAttribute('hidden', true);
        document.getElementById('gg_re_col').removeAttribute('hidden');
        document.getElementById('gg_lig_col').removeAttribute('hidden');
    }
});

document.getElementById('pcr_mmix_check').addEventListener('click', function()
{
    if(this.checked)
    {
        document.getElementById('pcr_mmix_col').removeAttribute('hidden');

    }
    else
    {
        document.getElementById('pcr_mmix_col').setAttribute('hidden', true);
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

// gibson
document.getElementById('gib_exo_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_gib_exonuclease_cost').removeAttribute('hidden');
        document.getElementById('id_gib_exonuclease_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_gib_exonuclease_cost').setAttribute('hidden', true);
        document.getElementById('id_gib_exonuclease_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('gib_ligase_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_gib_ligase_cost').removeAttribute('hidden');
        document.getElementById('id_gib_ligase_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_gib_ligase_cost').setAttribute('hidden', true);
        document.getElementById('id_gib_ligase_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('gib_poly_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_gib_polymerase_cost').removeAttribute('hidden');
        document.getElementById('id_gib_polymerase_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_gib_polymerase_cost').setAttribute('hidden', true);
        document.getElementById('id_gib_polymerase_n_reacts').setAttribute('hidden', true);
    }
});

// goldengate
document.getElementById('gg_re_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_gg_re_cost').removeAttribute('hidden');
        document.getElementById('id_gg_re_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_gg_re_cost').setAttribute('hidden', true);
        document.getElementById('id_gg_re_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('gg_ligase_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_gg_ligase_cost').removeAttribute('hidden');
        document.getElementById('id_gg_ligase_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_gg_ligase_cost').setAttribute('hidden', true);
        document.getElementById('id_gg_ligase_n_reacts').setAttribute('hidden', true);
    }
});

// // slic
// document.getElementById('slic_exo_check').addEventListener('click', function() {
//     if (this.checked) {
//         document.getElementById('id_slic_exonuclease_cost').removeAttribute('hidden');
//         document.getElementById('id_slic_exonuclease_n_reacts').removeAttribute('hidden');
//     }
//     else {
//         document.getElementById('id_slic_exonuclease_cost').setAttribute('hidden', true);
//         document.getElementById('id_slic_exonuclease_n_reacts').setAttribute('hidden', true);
//     }
// });

// document.getElementById('slic_ligase_check').addEventListener('click', function() {
//     if (this.checked) {
//         document.getElementById('id_slic_ligase_cost').removeAttribute('hidden');
//         document.getElementById('id_slic_ligase_n_reacts').removeAttribute('hidden');
//     }
//     else {
//         document.getElementById('id_slic_ligase_cost').setAttribute('hidden', true);
//         document.getElementById('id_slic_ligase_n_reacts').setAttribute('hidden', true);
//     }
// });

// biobricks
document.getElementById('bb_ligase_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_bb_ligase_cost').removeAttribute('hidden');
        document.getElementById('id_bb_ligase_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_bb_ligase_cost').setAttribute('hidden', true);
        document.getElementById('id_bb_ligase_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('bb_ecori_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_bb_EcoRI_cost').removeAttribute('hidden');
        document.getElementById('id_bb_EcoRI_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_bb_EcoRI_cost').setAttribute('hidden', true);
        document.getElementById('id_bb_EcoRI_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('bb_xbai_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_bb_XbaI_cost').removeAttribute('hidden');
        document.getElementById('id_bb_XbaI_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_bb_XbaI_cost').setAttribute('hidden', true);
        document.getElementById('id_bb_XbaI_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('bb_spei_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_bb_SpeI_cost').removeAttribute('hidden');
        document.getElementById('id_bb_SpeI_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_bb_SpeI_cost').setAttribute('hidden', true);
        document.getElementById('id_bb_SpeI_n_reacts').setAttribute('hidden', true);
    }
});

document.getElementById('bb_psti_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_bb_PstI_cost').removeAttribute('hidden');
        document.getElementById('id_bb_PstI_n_reacts').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_bb_PstI_cost').setAttribute('hidden', true);
        document.getElementById('id_bb_PstI_n_reacts').setAttribute('hidden', true);
    }
});