from django.shortcuts import render

def home(request):
    return render(request, 'assembly/home.html')

def about(request):
    return render(request, 'assembly/about.html', {'title': 'Assembly About'}) 

def gibson_form(request):
    return render(request, 'assembly/gibson_form.html', {'title': 'Gibson Assembly'})
    
def goldengate_form(request):
    return render(request, 'assembly/goldengate_form.html', {'title': 'Golden Gate Assembly'})