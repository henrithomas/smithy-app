# Generated by Django 3.2.5 on 2021-07-13 19:40

from django.db import migrations, models
import django.utils.timezone


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='GibsonAssembly',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(max_length=250)),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('backbone', models.CharField(max_length=10000)),
                ('insert', models.CharField(max_length=10000)),
                ('multipart', models.BooleanField(default=False)),
                ('addgene', models.BooleanField(default=False)),
                ('igem', models.BooleanField(default=False)),
                ('dnasu', models.BooleanField(default=False)),
                ('min_blast', models.PositiveIntegerField(verbose_name='min BLAST seq count')),
                ('max_blast', models.PositiveIntegerField(verbose_name='max BLAST seq count')),
                ('min_synth', models.PositiveIntegerField(verbose_name='min synthetic seq count')),
                ('max_synth', models.PositiveIntegerField(verbose_name='max synthetic seq count')),
                ('mv_conc', models.FloatField(verbose_name='monovalent ion concentration')),
                ('dv_conc', models.FloatField(verbose_name='divalent ion concentration')),
                ('dntp_conc', models.FloatField(verbose_name='dNTP concentration')),
                ('dna_conc', models.FloatField(verbose_name='dna concentration')),
                ('tm', models.FloatField(verbose_name='melting temperature')),
                ('overlap', models.PositiveIntegerField()),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='GoldenGateAssembly',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(max_length=250)),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('backbone', models.CharField(max_length=10000)),
                ('insert', models.CharField(max_length=10000)),
                ('multipart', models.BooleanField(default=False)),
                ('addgene', models.BooleanField(default=False)),
                ('igem', models.BooleanField(default=False)),
                ('dnasu', models.BooleanField(default=False)),
                ('min_blast', models.PositiveIntegerField(verbose_name='min BLAST seq count')),
                ('max_blast', models.PositiveIntegerField(verbose_name='max BLAST seq count')),
                ('min_synth', models.PositiveIntegerField(verbose_name='min synthetic seq count')),
                ('max_synth', models.PositiveIntegerField(verbose_name='max synthetic seq count')),
                ('mv_conc', models.FloatField(verbose_name='monovalent ion concentration')),
                ('dv_conc', models.FloatField(verbose_name='divalent ion concentration')),
                ('dntp_conc', models.FloatField(verbose_name='dNTP concentration')),
                ('dna_conc', models.FloatField(verbose_name='dna concentration')),
                ('tm', models.FloatField(verbose_name='melting temperature')),
                ('overhangs', models.IntegerField(choices=[(15, '15-overhangs'), (20, '20-overhangs'), (25, '25-overhangs'), (30, '30-overhangs')], verbose_name='overhang count')),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
