module Mandelbrot.Graphics.Text where

import Control.Exception (finally)
import Data.Complex (Complex((:+)))
import Data.Word (Word8)
import Numeric (showEFloat)

import Mandelbrot.Graphics.Colour
import Mandelbrot.Graphics.Compute
import Mandelbrot.Graphics.Render
import Mandelbrot.Graphics.Transform

ascii1 :: Compute Double -> Char
ascii1 Unknown{} = '?'
ascii1 e = "@#*+;:,. " !! round (0 `max` (7 + logBase' 2 (computeDE e)) `min` 8)
  where
    logBase' b x
      | x <= 0 = -1 / 0
      | otherwise = logBase b x

ascii :: [[Compute Double]] -> String
ascii = unlines . map (map ascii1) . map double
  where
    double (a:as) = a : a : double as
    double _ = []

type RGB8 = (Word8, Word8, Word8)

blend8 :: RGB8 -> RGB8 -> Double -> RGB8
blend8 (r1,g1,b1) (r2,g2,b2) t = (f r1 r2, f g1 g2, f b1 b2)
  where
   f x y = round $ 255 * sqrt (blend ((fromIntegral x / 255)^(2::Int)) ((fromIntegral y / 255)^(2::Int)) t)

blend :: Double -> Double -> Double -> Double
blend a b t = a + t * (b - a)

fg, bg :: RGB8 -> String
fg (r,g,b) = "\ESC[38;2;" ++ show r ++ ";" ++ show g ++ ";" ++ show b ++ "m"
bg (r,g,b) = "\ESC[48;2;" ++ show r ++ ";" ++ show g ++ ";" ++ show b ++ "m"

ansi1 :: Bool -> (Word8, Word8, Word8) -> (Word8, Word8, Word8) -> String
ansi1 True f b = fg f ++ bg b ++ [upperHalfBlock]
ansi1 False f b = bg (blend8 f b 0.5) ++ [' ']

upperHalfBlock :: Char
upperHalfBlock = '\x2580'

ansi :: Bool -> [[(Word8, Word8, Word8)]] -> String
ansi u (a:b:cs) = concat (zipWith (ansi1 u) a b) ++ "\ESC[m\n" ++ ansi u cs
ansi u [a] = concat (zipWith (ansi1 u) a (repeat (128,128,128))) ++ "\ESC[m\n"
ansi _ [] = ""

putImage_ :: Bool -> Complex Double -> Double -> Int -> IO ()
putImage_ u c@(x:+y) r n = flip finally reset $ do
  putStr . ansi u . map (map (grayscale . monochrome)) $
    scanline w h (rectangular w h c r) 1e6 n
  putStrLn $ show x ++ " + " ++ show y ++ " i @ " ++
               showEFloat (Just 3) r "" ++
               " (" ++ show n ++ " iterations)"
  where
    w = 80
    h = 38

putImage', putImage :: Complex Double -> Double -> Int -> IO ()
putImage' = putImage_ False
putImage = putImage_ True

reset :: IO ()
reset = putStr "\ESC[m"
