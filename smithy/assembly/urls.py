from collections import namedtuple
from django.urls import path
from . import views
from .views import (
    GibsonDetailView,
    GibsonCreateView,
    GibsonSolutionDetailView,
    GibsonPartDetailView,
    GibsonPrimerDetailView,
    GoldenGateDetailView,
    GoldenGateCreateView,
    GoldenGateSolutionDetailView,    
    GoldenGatePartDetailView,
    GoldenGatePrimerDetailView,
)

urlpatterns = [
    path('', views.home, name='assembly-home'),
    path('about/', views.about, name='assembly-about'),
    path('gibson/<int:pk>/', GibsonDetailView.as_view(), name='gibson-detail'),
    path('gibson/new/', GibsonCreateView.as_view(), name='gibson-create'),
    path('gibson/solution/<int:pk>', GibsonSolutionDetailView.as_view(), name='gibson-solution-detail'),
    path('gibson/part/<int:pk>', GibsonPartDetailView.as_view(), name='gibson-part-detail'),
    path('gibson/primer/<int:pk>', GibsonPrimerDetailView.as_view(), name='gibson-primer-detail'),
    path('goldengate/<int:pk>/', GoldenGateDetailView.as_view(), name='goldengate-detail'),
    path('goldengate/new/', GoldenGateCreateView.as_view(), name='goldengate-create'),
    path('goldengate/solution/<int:pk>', GoldenGateSolutionDetailView.as_view(), name='goldengate-solution-detail'),
    path('goldengate/part/<int:pk>', GoldenGatePartDetailView.as_view(), name='goldengate-part-detail'),
    path('goldengate/primer/<int:pk>', GoldenGatePrimerDetailView.as_view(), name='goldengate-primer-detail'),
]