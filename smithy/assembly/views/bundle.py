from ..models import AssemblyBundle
from django.views.generic import (
    DetailView,
    ListView
)
from django.shortcuts import render, redirect
from ..forms import BundleForm
from ..services.bundle import bundle_create_service

class AssemblyBundleDetailView(DetailView):
    model = AssemblyBundle
    context_object_name = 'assembly_bundle'

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)
        # TODO change how the queries happen here, maybe change to code in the GET request
        context['title'] = self.object.title 
        
        context['gibson'] = self.object.gibson.filter().first()
        if context['gibson']:
            context['gibson_solutions'] = context['gibson'].gibsonsolution_set.all()
        
        context['goldengate'] = self.object.goldengate.filter().first()
        if context['goldengate']:
            context['goldengate_solutions'] = context['goldengate'].goldengatesolution_set.all()

        context['slic'] = self.object.slic.filter().first()
        if context['slic']:
            context['slic_solutions'] = context['slic'].slicsolution_set.all()
        
        context['pcr'] = self.object.pcr.filter().first()
        if context['pcr']:
            context['pcr_solutions'] = context['pcr'].pcrsolution_set.all()

        context['biobrick'] = self.object.biobricks.filter().first()
        if context['biobrick']:
            context['biobrick_solutions'] = context['biobrick'].biobrickssolution_set.all()

        return context


class AssemblyBundleListView(ListView):
    model = AssemblyBundle
    context_object_name = 'bundles'

    def get_queryset(self):
        qs = super().get_queryset()
        return qs.all() \
            .only('title', 'date_created', 'description') \
            .order_by('-date_created')
    
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'Assembly Bundles'
        return context


def assembly_bundle(request):
    if request.method == 'POST':
        bundle_form = BundleForm(request.POST, request.FILES)

        if bundle_form.is_valid():
            bundle_pk = bundle_create_service(bundle_form.cleaned_data)
            return redirect('bundle-detail', bundle_pk)
    else:
        bundle_form = BundleForm()
    return render(request, 'assembly/assemblybundle_form.html', {'bundle_form': bundle_form})