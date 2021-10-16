from django import template

register = template.Library()

@register.simple_tag
def assembly_obj_url(assembly):
    if assembly.__class__.__name__ == 'GibsonAssembly':
        return f'https://smithy.fungalgenomics.ca/assembly/gibson/{assembly.pk}/'
    elif assembly.__class__.__name__ == 'GoldenGateAssembly':
        return f'https://smithy.fungalgenomics.ca/assembly/goldengate/{assembly.pk}/'
    elif assembly.__class__.__name__ == 'BioBricksAssembly':
        return f'https://smithy.fungalgenomics.ca/assembly/biobricks/{assembly.pk}/'
    elif assembly.__class__.__name__ == 'PCRAssembly':
        return f'https://smithy.fungalgenomics.ca/assembly/pcr/{assembly.pk}/'
    else:
        return f'https://smithy.fungalgenomics.ca/assembly/slic/{assembly.pk}/'

@register.simple_tag
def assembly_type(assembly):
    if assembly.__class__.__name__ == 'GibsonAssembly':
        return 'Gibson'
    elif assembly.__class__.__name__ == 'GoldenGateAssembly':
        return 'Golden Gate'
    elif assembly.__class__.__name__ == 'BioBricksAssembly':
        return 'BioBricks'
    elif assembly.__class__.__name__ == 'PCRAssembly':
        return 'PCR'
    else:
        return f'SLIC'