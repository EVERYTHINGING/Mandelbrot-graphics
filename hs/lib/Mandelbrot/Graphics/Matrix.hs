module Mandelbrot.Graphics.Matrix
  ( Mat2(..)
  , identity
  , det
  , tr
  , inv
  , mul
  , diagonalize
  , moebius3
  , interpolate
  ) where

data Mat2 r = Mat2 !r !r !r !r

identity :: Num r => Mat2 r
identity = Mat2 1 0 0 1

det, tr :: Num r => Mat2 r -> r
det (Mat2 a b c d) = a * d - b * c
tr (Mat2 a _ _ d) = a + d

inv :: Fractional r => Mat2 r -> Mat2 r
inv m@(Mat2 a b c d) =
  let e = det m
  in  Mat2 (d / e) (-b/e) (-c/e) (a/e)

mul :: Num r => Mat2 r -> Mat2 r -> Mat2 r
mul (Mat2 la lb lc ld) (Mat2 ra rb rc rd) = Mat2 a b c d
  where
    a = la * ra + lb * rc
    b = la * rb + lb * rd
    c = lc * ra + ld * rc
    d = lc * rb + ld * rd

diagonalize :: (Eq r, Floating r) => Mat2 r -> (Mat2 r, Mat2 r, Mat2 r)
diagonalize m@(Mat2 ma mb mc md) = (p, d, p1)
  where
    d = Mat2 l1 0 0 l2
    l1 = tr2 + k
    l2 = tr2 - k
    k = sqrt (tr2 * tr2 - det m)
    tr2 = tr m / 2
    p1 = inv p
    p | mb /= 0 = Mat2 mb mb (l1 - ma) (l2 - ma)
      | mc /= 0 = Mat2 (l1 - md) (l2 - md) mc mc
      | otherwise = identity

moebius3 :: Num r => r -> r -> r -> Mat2 r
moebius3 zero one infinity = Mat2 a b c d
  where
    a = infinity * (zero - one)
    b = zero * (one - infinity)
    c = zero - one
    d = one - infinity

interpolate :: (Eq r, Floating r) => Mat2 r -> Mat2 r -> r -> Mat2 r
interpolate f g = \t ->
    let e = Mat2 (da ** t) 0 0 (dd ** t)
        fpe = mul fp e
    in  mul fpe p1
  where
    f1 = inv f
    f1g = mul f1 g
    (p, Mat2 da _ _ dd, p1) = diagonalize f1g
    fp = mul f p
