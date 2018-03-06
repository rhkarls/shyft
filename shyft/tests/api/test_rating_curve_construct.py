import unittest
from shyft.api import RatingCurveSegment, RatingCurveSegments, RatingCurveFunction, RatingCurveTimeFunction, RatingCurveTimeFunctions, RatingCurveParameters, Calendar


class FlowRatingCurves(unittest.TestCase):

    def test_easy_rating_curve_construct(self):
        utc = Calendar()
        rating_curve = RatingCurveParameters(RatingCurveTimeFunctions([
            RatingCurveTimeFunction(
                utc.time(1950, 3, 27), RatingCurveFunction(RatingCurveSegments([
                    RatingCurveSegment(lower=0.474, a=5.97489, b=-0.4745, c=2.36997)
                ]))),
            RatingCurveTimeFunction(
                utc.time(1968, 7, 29), RatingCurveFunction(RatingCurveSegments([
                    RatingCurveSegment(lower=0.25, a=2.9822, b=-0.45, c=1.5078),
                    RatingCurveSegment(lower=0.79, a=3.9513, b=-0.45, c=2.8087),
                    RatingCurveSegment(lower=1.38, a=5.7071, b=-0.45, c=2.3503),
                    RatingCurveSegment(lower=2.55, a=8.2672, b=-0.45, c=2.052)
                ])))
        ]))
        self.assertIsNotNone(rating_curve)
        flow = rating_curve.flow(utc.time(2018, 1, 1), 3.2)
        self.assertAlmostEqual(flow, 117.8103380205204)
