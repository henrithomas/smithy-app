from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='assembly-home'),
    path('about/', views.about, name='assembly-about'),
    path('gibson/', views.gibson_form, name='gibson-form'),
    path('goldengate/', views.goldengate_form, name='goldengate-form'),
    path('goldengate_dev/', views.goldengate_dev, name='goldengate-dev'),
    path('assembly_submit/', views.submission, name='assembly-submit'),
]