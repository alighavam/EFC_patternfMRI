import os
import numpy
from cfmm_base import infotodict as cfmminfodict

def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return (template, outtype, annotation_classes)

def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    allowed template fields - follow python string module:
    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    
    # call cfmm for general labelling and get dictionary
    # info = cfmminfodict(seqinfo)
    bold_magnitude = create_key('{bids_subject_session_dir}/func/{bids_subject_session_prefix}_run-{item:02d}_bold_magnitude')
    bold_phase = create_key('{bids_subject_session_dir}/func/{bids_subject_session_prefix}_run-{item:02d}_bold_phase')
    sbref = create_key('{bids_subject_session_dir}/func/{bids_subject_session_prefix}_run-{item:02d}_sbref')
    
    # fmap_sbref = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_dir-{dir}_epi')

    #GRE phase diff 
    fmap_diff = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_phasediff')
    fmap_magnitude = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_magnitude')


    # t13d = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-MPRAGE_run-{item:02d}_T1w')
    info = {bold_phase:[],
            bold_magnitude:[],
            sbref:[],
            fmap_diff:[],
            fmap_magnitude:[]}

    for idx, s in enumerate(seqinfo):
        if ('bold' in (s.series_description).strip()):
            if ('SBRef' in (s.series_description).strip()):
                info[sbref].append({'item': s.series_id})
            elif (s.dim4 > 300):
                if('P' in (s.image_type[2].strip()) ):
                    info[bold_phase].append({'item': s.series_id})
                if('M' in (s.image_type[2].strip()) ):
                    info[bold_magnitude].append({'item': s.series_id})

        if ('field_mapping' in s.protocol_name):   
            if (s.dim4 == 1) and ('gre_field_mapping' == (s.series_description).strip()):
                if('P' in (s.image_type[2].strip()) ):
                    info[fmap_diff].append({'item': s.series_id})
                if('M' in (s.image_type[2].strip()) ):
                    info[fmap_magnitude].append({'item': s.series_id})

        # elif ('1.8iso_PA' in (s.series_description).strip()):
        #    if (s.dim4 == 1):
                # if 'SBRef' in (s.series_description).strip():
        #        info[fmap_sbref].append({'item': s.series_id, 'dir': 'PA'})

        # elif ('mp2rage' in (s.series_description).strip()):
        #    info[t13d].append({'item': s.series_id})

    return info