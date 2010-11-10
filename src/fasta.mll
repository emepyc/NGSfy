(************************************************************************)
(*
   This file is part of NGSfy software. NGSfy is free software; you can
   redistribute it and/or modify it under the terms of the GNU General
   Public License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*)
(***********************************************************************)
(* fasta.mll : Input of fasta sequences *)

{
open ExtLib
}

rule pfasta = parse
| '#' ([^'\r''\n'])* ['\r''\n'] { pfasta lexbuf }
   |'>'([^'\r' '\n'])* as header {
      let sequence = getsequence lexbuf in ((String.lchop header),(MyString.cleanseq sequence)) }
   | eof { raise End_of_file }
 and getsequence = parse
   |([^'>'])* as seq { seq }


{}
