from django.shortcuts import render,redirect
from django.views.generic.edit import CreateView
from .forms import GoldenGateForm, GibsonForm
from .models import GibsonAssembly, GoldenGateAssembly
from django.views.generic import (
    DetailView,
    CreateView
)
from django.contrib.messages.views import SuccessMessageMixin

from assemblies.gibson import GibsonAssembler

def db_list(addgene, igem, dnasu):
    l = []
    if addgene:
        l.append('addgene')
    if igem: 
        l.append('igem')
    if dnasu:
        l.append('dnasu')
    return l

def home(request):
    return render(request, 'assembly/home.html')

def about(request):
    return render(request, 'assembly/about.html', {'title': 'Assembly About'}) 

def submission(request):
    return render(request, 'assembly/submission.html', {'title': 'Assembly Submission'}) 

def gibson_form(request):
    return render(request, 'assembly/gibson_form.html', {'title': 'Gibson Assembly'})
    
def goldengate_form(request):
    return render(request, 'assembly/goldengate_form.html', {'title': 'Golden Gate Assembly'})

def goldengate_dev(request):
    if request.method == 'POST':
        form = GoldenGateForm(request.POST)
        if form.is_valid():
            return redirect('assembly-submit') 
    else:
        form = GoldenGateForm()

    return render(request, 'assembly/goldengate_dev.html', {'title': 'Golden Gate Assembly', 'form': form})

def gibson_dev(request):
    if request.method == 'POST':
        form = GibsonForm(request.POST)
        if form.is_valid():
            return redirect('assembly-submit') 
    else:
        form = GibsonForm()

    return render(request, 'assembly/gibson_dev.html', {'title': 'Gibson Assembly', 'form': form})

class GibsonDetailView(DetailView):
    model = GibsonAssembly
    context_object_name = 'gibson'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.title
        return context

class GibsonCreateView(SuccessMessageMixin, CreateView):
    model = GibsonAssembly
    success_message = 'View your new Gibson assembly below...'
    fields = [
            'title',
            'backbone_file',
            'insert_file',
            'addgene',
            'igem',
            'dnasu',
            'min_blast',
            'max_blast',
            'min_synth',
            'max_synth',
            'mv_conc',
            'dv_conc',
            'dntp_conc',
            'dntp_conc',
            'dna_conc', 
            'tm',
            'overlap']

    def form_valid(self, form):
        self.object = form.save()
        gib_assembler = GibsonAssembler(
                            self.object.mv_conc, 
                            self.object.dv_conc, 
                            self.object.dna_conc,
                            self.object.dntp_conc, 
                            self.object.tm, 
                            self.object.backbone_file.path, 
                            self.object.insert_file.path, 
                            db_list(self.object.addgene, self.object.igem, self.object.dnasu), 
                            min_frag=self.object.min_blast, 
                            max_frag=self.object.max_blast, 
                            min_synth=self.object.min_synth, 
                            max_synth=self.object.max_synth)
        # results, error = gib_assembler.query()
        # gib_assembler.solution_building(results)
        # gib_assembly, gib_fragments = gib_assembler.design(solution=1)
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'New Gibson Assembly'
        return context


class GoldenGateDetailView(DetailView):
    model = GoldenGateAssembly
    context_object_name = 'goldengate'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.title
        return context

class GoldenGateCreateView(SuccessMessageMixin, CreateView):
    model = GoldenGateAssembly
    success_message = 'View your new Golden Gate assembly below...'
    fields = [
            'title',
            'backbone_file',
            'insert_file',
            'addgene',
            'igem',
            'dnasu',
            'min_blast',
            'max_blast',
            'min_synth',
            'max_synth',
            'mv_conc',
            'dv_conc',
            'dntp_conc',
            'dntp_conc',
            'dna_conc', 
            'tm',
            'overhangs']

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'New Golden Gate Assembly'
        return context