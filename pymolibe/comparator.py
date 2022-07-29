from itertools import combinations
from numpy import ones, dot, transpose, linalg, sqrt, power, sum, abs, max


def align(candidate_structure, reference_structure):
    """
    Rotate candidate structure unto reference structure using Kabsch algorithm based on each position.

    :param candidate_structure: candidate structure represented by three-dimension position list.
    :type candidate_structure: numpy.ndarray

    :param reference_structure: reference structure represented by three-dimension position list.
    :type reference_structure: numpy.ndarray

    :returns candidate structure, reference structure: final candidate structure and reference structure.
    :rtype numpy.ndarray, numpy.ndarray

    .. note::
        reference [1]: Yang Zhang and Jeffrey Skolnick (2004) Proteins

        reference [2]: Wolfgang Kabsch (1976) Acta Crystallogr. D.
    """
    for center_index in range(len(candidate_structure)):
        # unify the center point.
        candidate = candidate_structure - candidate_structure[center_index]
        reference = reference_structure - reference_structure[center_index]
        center = dot(transpose(candidate), reference)

        # decompose the singular value for rotation matrix.
        translation, _, unitary = linalg.svd(center)
        if (linalg.det(translation) * linalg.det(unitary)) < 0.0:
            translation[:, -1] = -translation[:, -1]

        # create rotation matrix and rotate source structure.
        rotation = dot(translation, unitary)
        candidate = dot(candidate, rotation)

        # recover the original state of reference structure.
        candidate += reference_structure[center_index]
        reference += reference_structure[center_index]

        yield candidate, reference


class Similarity(object):

    def __init__(self, name, threshold=0.1, is_state=True):
        """
        Initialize a similarity method.

        :param name: name of the similarity calculation.
        :type name: str

        :param threshold: threshold of the similarity evaluation.
        :type threshold: float

        :param is_state: output is similarity of not; otherwise, output is the score or value of similarity.
        :type is_state: bool

        .. note::
            The default threshold is 0.1, from GAxel T. Brünger (1992) Nature
        """
        self.name = name
        self.threshold = threshold
        self.is_state = is_state

    def calculate(self, structure_1, resolution_1, structure_2, resolution_2):
        """
        Calculate the similarity between two structures.

        :param structure_1: one structure represented by three-dimension position list one.
        :type structure_1: numpy.ndarray

        :param resolution_1: one resolution of position list one.
        :type resolution_1: float

        :param structure_2: another structure represented by three-dimension position list two.
        :type structure_2: numpy.ndarray

        :param resolution_2: another resolution of position list two.
        :type resolution_2: float

        :return: If the score more than or the difference is less than threshold function, return True, otherwise False.
        :rtype dict
        """
        raise NotImplementedError

    def __call__(self, structure_1, resolution_1, structure_2, resolution_2):
        return self.calculate(structure_1, resolution_1, structure_2, resolution_2)


class GMD(Similarity):

    def __init__(self, threshold=0.1, is_state=True):
        """
        Calculate similarity by global maximum distance.

        :param threshold: maximum difference threshold.
        :type threshold: float

        :param is_state: output is similarity of not; otherwise, output is the score or value of similarity.
        :type is_state: bool
        """
        super().__init__(name="GMD", threshold=threshold, is_state=is_state)

    def calculate(self, structure_1, resolution_1, structure_2, resolution_2):
        if structure_1.shape != structure_2.shape:
            if self.is_state:
                return False
            else:
                return -ones(shape=(len(structure_1) if len(structure_1) > len(structure_2) else len(structure_2),))

        error_range, values = (resolution_1 + resolution_2) * self.threshold, []

        # calculate the distance between any two corresponding vertices.
        for index_1, index_2 in combinations(list(range(len(structure_1))), 2):
            length_1 = linalg.norm(structure_1[index_1] - structure_1[index_2], ord=2)
            length_2 = linalg.norm(structure_2[index_1] - structure_2[index_2], ord=2)
            values.append(abs(length_1 - length_2))
            if self.is_state and abs(length_1 - length_2) > error_range:
                return False

        return max(values) > error_range if self.is_state else values


class GMV(Similarity):

    def __init__(self, threshold=0.1, is_state=True):
        """
        Calculate similarity by global maximum vector.

        :param threshold: maximum difference threshold.
        :type threshold: float

        :param is_state: output is similarity of not; otherwise, output is the score or value of similarity.
        :type is_state: bool
        """
        super().__init__(name="GMV", threshold=threshold, is_state=is_state)

    def calculate(self, structure_1, resolution_1, structure_2, resolution_2):
        if structure_1.shape != structure_2.shape:
            if self.is_state:
                return False
            else:
                return -ones(shape=(len(structure_1) if len(structure_1) > len(structure_2) else len(structure_2),))

        error_range, values = (resolution_1 + resolution_2) * self.threshold, []

        # calculate any two corresponding vector with length 3.
        for index in range(len(structure_1) - 4 + 1):
            candidate, reference = align(structure_1[index: index + 4], structure_2[index: index + 4])
            value = sqrt(sum(power(candidate - reference, 2)) / len(structure_1))
            values.append(value)
            if self.is_state and value > error_range:
                return False

        return max(values) > error_range if self.is_state else values


class RMSD(Similarity):

    def __init__(self, threshold=0.1, is_state=True):
        """
        Calculate similarity by root-mean-square deviation.

        :param threshold: maximum difference threshold.
        :type threshold: float

        :param is_state: output is similarity of not; otherwise, output is the score or value of similarity.
        :type is_state: bool

        .. note::
            reference: Evangelos A Coutsias, Chaok Seok, and Ken A Dill (2004) J. Comput. Chem.
        """
        super().__init__(name="RMSD", threshold=threshold, is_state=is_state)

    def calculate(self, structure_1, resolution_1, structure_2, resolution_2):
        if structure_1.shape != structure_2.shape:
            return False if self.is_state else -1

        error_range, minimum_distance, structures = (resolution_1 + resolution_2) * self.threshold, None, None

        for candidate, reference in align(structure_1, structure_2):
            # calculate the value of root-mean-square deviation.
            distance = sqrt(sum(power(candidate - reference, 2)) / len(structure_1))
            if minimum_distance is None or minimum_distance > distance:
                structures, minimum_distance = (candidate, reference), distance

        return minimum_distance <= error_range if self.is_state else minimum_distance


class TMScore(Similarity):

    def __init__(self, threshold=0.7, is_state=True):
        """
        Calculate similarity by the template modeling score (TM-score).

        :param threshold: score minimum threshold.
        :type threshold: float

        :param is_state: output is similarity of not; otherwise, output is the score or value of similarity.
        :type is_state: bool

        .. note::
            Score reference: Yang Zhang and Jeffrey Skolnick (2004) Proteins

            The default threshold is 0.5, from Jinrui Xu and Yang Zhang (2010) Bioinformatics
            The default threshold is 0.7, from Andrew W Senior, Richard Evans, John Jumper, and others (2020) Nature

            Comments from Yang Zhang's lab:
            In our study ‘How significant is a protein structure similarity with TM-score = 0.5?’,
            we only obtain statistics for structures and fragments with more than 20 residues.
            Whether or not fragments shorter than 20 residues,
            which is the focus of your study, follow the same statistics is unknown.
        """
        super().__init__(name="TM-score", threshold=threshold, is_state=is_state)

    def calculate(self, structure_1, resolution_1, structure_2, resolution_2):
        if structure_1.shape != structure_2.shape:
            return False if self.is_state else 0.0

        # calculate the distance scale that normalizes distances.
        distance_scale = 1.24 * (len(structure_1) - 15) ** (1.0 / 3.0) - 1.8

        # rotate source structure unto target structure.
        maximum_score, structures = 0, None
        for candidate, reference in align(structure_1, structure_2):
            distances = linalg.norm(candidate - reference, ord=2, axis=1)
            score = sum(1.0 / (1.0 + (distances / distance_scale) ** 2)) / len(structure_1)
            if maximum_score < score:
                structures, maximum_score = (candidate, reference), score

        return maximum_score >= self.threshold if self.is_state else maximum_score


class GDTScore(Similarity):

    def __init__(self, threshold=90.0, is_state=True, is_total_score=True):
        """
        Calculate similarity by the global distance test score (GDT-score).

        :param threshold: score minimum threshold.
        :type threshold: float.

        :param is_state: output is similarity of not; otherwise, output is the score or value of similarity.
        :type is_state: bool

        :param is_total_score: score type is total score (GDT_TS); otherwise, score type is high accurate (GDT_HA).
        :type is_total_score: bool.

        .. note::
            Score reference: Adam Zemla (2003) Nucleic Acids Res.
            The default threshold is 90.0, from Gaurav Chopra, Nir Kalisman, and Michael Levitt (1997) Proteins
        """
        super().__init__(name="GDT-score", threshold=threshold, is_state=is_state)
        self.is_total_score = is_total_score

    def calculate(self, structure_1, resolution_1, structure_2, resolution_2):
        if structure_1.shape != structure_2.shape:
            return False if self.is_state else 0.0

        scores, length = [], len(structure_1)
        for candidate, reference in align(structure_1, structure_2):
            # calculate the value of root-mean-square deviation.
            values = sqrt(sum(power(candidate - reference, 2), axis=1))

            # count the distances under different selected cut offs.
            if self.is_total_score:
                cutoffs = {"1": 0, "2": 0, "4": 0, "8": 0}
            else:
                cutoffs = {"0.5": 0, "1": 0, "2": 0, "4": 0}

            for value in values:
                for cutoff in cutoffs.keys():
                    if value <= float(cutoff):
                        cutoffs[cutoff] += 1

            # calculate the final score.
            scores.append(sum(cutoffs.values()) / (4.0 * length) * 100.0)

        return max(scores) >= self.threshold if self.is_state else max(scores)
