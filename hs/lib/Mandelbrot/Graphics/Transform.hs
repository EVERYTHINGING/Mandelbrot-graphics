module Mandelbrot.Graphics.Transform
  ( Transform(..)
  , invert
  , compose
  , rectangular
  , linear
  , moebius
  , moebius3
  , moebius6
  , cardioid
  ) where

import Prelude hiding (reverse)

import Data.Complex (Complex((:+)), conjugate)

import Mandelbrot.Graphics.Matrix (Mat2(..))
import qualified Mandelbrot.Graphics.Matrix as M

data Transform r = Transform
  { forward, reverse :: (Complex r, Complex r) -> (Complex r, Complex r) }

invert :: Transform r -> Transform r
invert (Transform f r) = Transform r f

compose :: Transform r -> Transform r -> Transform r
compose (Transform a a') (Transform b b') = Transform (b . a) (a' . b')

rectangular :: RealFloat r => Int -> Int -> Complex r -> r -> Transform r
rectangular width height c r = Transform
  { forward = \(c0, dc0) ->
    ( c + ((r / h2) :+ 0) * conjugate (c0 - ((w2 - 0.5) :+ (h2 - 0.5)))
    , conjugate ((r / (fromIntegral height / 2) :+ 0) * dc0)
    )
  , reverse = \(c0, dc0) ->
    ( conjugate ((c0 - c) * ((h2 / r) :+ 0)) + ((w2 - 0.5) :+ (h2 - 0.5))
    , conjugate (((h2 / r) :+ 0) * dc0)
    )
  }
  where
    h2 = fromIntegral height / 2
    w2 = fromIntegral width / 2

linear :: RealFloat r => Complex r -> Complex r -> Transform r
linear mul add = Transform
  { forward = \(c0, dc0) -> (mul * c0 + add, mul * dc0)
  , reverse = \(c0, dc0) -> ((c0 - add) / mul, dc0 / mul)
  }

moebius :: RealFloat r => Mat2 (Complex r) -> Transform r
moebius m@(Mat2 a b c d) = Transform
  { forward = \(c0, dc0) ->
      let e = c * c0 + d
      in  ((a * c0 + b) / e, dc0 * M.det m / (e * e))
  , reverse = \(c0, dc0) ->
      let e = a - c * c0
      in  ((d * c0 - b) / e, dc0 * M.det m / (e * e))
  }

moebius3 :: RealFloat r => Complex r -> Complex r -> Complex r -> Transform r
moebius3 zero one infinity = moebius (M.moebius3 zero one infinity)

moebius6 :: RealFloat r => Complex r -> Complex r -> Complex r -> Complex r -> Complex r -> Complex r -> Transform r
moebius6 a0 a1 a2 b0 b1 b2 = moebius (M.mul (M.moebius3 a0 a1 a2) (M.inv (M.moebius3 b0 b1 b2)))

cardioid :: RealFloat r => Transform r
cardioid = Transform
  { forward = \(c0, dc0) ->
      let r = sqrt (1 - 4 * c0)
      in  (c0 - 1, dc0 * (-2) / r)
  , reverse = \(c0, dc0) ->
      let r = c0 + 1
      in  ((1 - r * r) / 4, dc0 * r / (-2))
  }
