from django.shortcuts import render,redirect
from .forms import GoldenGateForm

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
            print(form.cleaned_data['addgene'])
            return redirect('assembly-submit') 
    else:
        form = GoldenGateForm()

    return render(request, 'assembly/goldengate_dev.html', {'title': 'Golden Gate Assembly', 'form': form})