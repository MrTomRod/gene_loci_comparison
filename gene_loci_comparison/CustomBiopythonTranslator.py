from copy import copy
from Bio.SeqFeature import SeqFeature
from dna_features_viewer import GraphicFeature, BiopythonTranslator


class CustomBiopythonTranslator(BiopythonTranslator):
    label_fields = None

    def __init__(self, label_fields: [str], *args, **kwargs):
        self.label_fields = label_fields
        super().__init__(*args, **kwargs)

    def compute_filtered_features(self, features):
        """Return the list of features minus the ignored ones.

        If a feature has multiple locations, it will be split into multiple
        features with a location.
        """
        filtered_features = []
        for f in features:
            if all([fl(f) for fl in self.features_filters]) and f.type not in self.ignored_features_types:
                if len(f.location.parts) == 1:
                    # simple case: only one location
                    filtered_features.append(f)
                else:
                    # multiple locations -> create multiple features
                    for loc in f.location.parts:
                        f_copy = copy(f)
                        f_copy.location = loc
                        filtered_features.append(f_copy)
        return filtered_features
