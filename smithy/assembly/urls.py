from django.urls import path
from .views.base import *
from .views.gibson import *
from .views.goldengate import *
from .views.pcr import *
from .views.slic import *
from .views.biobricks import *
from .views.bundle import *

urlpatterns = [
    path('', home, name='assembly-home'),
    path('all/', assemblies_list, name='assemblies-list'),
    path('about/', about, name='assembly-about'),
    path('bundle/', assembly_bundle, name='assembly-bundle'),
    path('bundle/<int:pk>/', AssemblyBundleDetailView.as_view(), name='bundle-detail'),
    path('bundle/list/', AssemblyBundleListView.as_view(), name='bundle-list'),
    path('gibson/<int:pk>/', GibsonDetailView.as_view(), name='gibson-detail'),
    path('gibson/new/', GibsonCreateView.as_view(), name='gibson-create'),
    path('gibson/list/', GibsonListView.as_view(), name='gibson-list'),
    path('gibson/solution/<int:pk>/', GibsonSolutionDetailView.as_view(), name='gibson-solution-detail'),
    path('gibson/part/<int:pk>/', GibsonPartDetailView.as_view(), name='gibson-part-detail'),
    path('gibson/primer/<int:pk>/', GibsonPrimerDetailView.as_view(), name='gibson-primer-detail'),
    path('goldengate/<int:pk>/', GoldenGateDetailView.as_view(), name='goldengate-detail'),
    path('goldengate/new/', GoldenGateCreateView.as_view(), name='goldengate-create'),
    path('goldengate/list/', GoldenGateListView.as_view(), name='goldengate-list'),
    path('goldengate/solution/<int:pk>/', GoldenGateSolutionDetailView.as_view(), name='goldengate-solution-detail'),
    path('goldengate/part/<int:pk>/', GoldenGatePartDetailView.as_view(), name='goldengate-part-detail'),
    path('goldengate/primer/<int:pk>/', GoldenGatePrimerDetailView.as_view(), name='goldengate-primer-detail'),
    path('biobricks/<int:pk>/', BioBricksDetailView.as_view(), name='biobricks-detail'),
    path('biobricks/new/', BioBricksCreateView.as_view(), name='biobricks-create'),
    path('biobricks/list/', BioBricksListView.as_view(), name='biobricks-list'),
    path('biobricks/solution/<int:pk>/', BioBricksSolutionDetailView.as_view(), name='biobricks-solution-detail'),
    path('biobricks/part/<int:pk>/', BioBricksPartDetailView.as_view(), name='biobricks-part-detail'),
    path('biobricks/primer/<int:pk>/', BioBricksPrimerDetailView.as_view(), name='biobricks-primer-detail'),
    path('pcr/<int:pk>/', PCRDetailView.as_view(), name='pcr-detail'),
    path('pcr/new/', PCRCreateView.as_view(), name='pcr-create'),
    path('pcr/list/', PCRListView.as_view(), name='pcr-list'),
    path('pcr/solution/<int:pk>/', PCRSolutionDetailView.as_view(), name='pcr-solution-detail'),
    path('pcr/part/<int:pk>/', PCRPartDetailView.as_view(), name='pcr-part-detail'),
    path('pcr/primer/<int:pk>/', PCRPrimerDetailView.as_view(), name='pcr-primer-detail'),
    path('slic/<int:pk>/', SLICDetailView.as_view(), name='slic-detail'),
    path('slic/new/', SLICCreateView.as_view(), name='slic-create'),
    path('slic/list/', SLICListView.as_view(), name='slic-list'),
    path('slic/solution/<int:pk>/', SLICSolutionDetailView.as_view(), name='slic-solution-detail'),
    path('slic/part/<int:pk>/', SLICPartDetailView.as_view(), name='slic-part-detail'),
    path('slic/primer/<int:pk>/', SLICPrimerDetailView.as_view(), name='slic-primer-detail'),
]
