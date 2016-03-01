{-# LANGUAGE DeriveDataTypeable #-}

{-| Define common data types for features of sequences -}

module Bio.Core.Strand (Strand(..)) where

import Data.Typeable (Typeable)

-- | A 'Strand' is either plus (forward) or minus (reverse or reverse-complement)
data Strand = Plus | Minus deriving (Eq,Ord,Read,Show,Typeable)
