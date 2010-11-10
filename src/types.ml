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
open Bin_io

let (~!) = (lnot);;
let (&!) = (land);;

let sff_magic = Int32.of_int 0x2e736666 (* ".sff" *)
let sff_version = "\000\000\000\001"

type sff_common_header = {
    magic:int32;
    version:string;
    index_offset:int64;
    index_len:int32;
    nreads:int32;
    header_len:int;
    key_len:int;
    flow_len:int;
    flowgram_format:int;
    mutable flow:string;
    mutable key:string
  } 

(*
 * We have one read_header per "reading" in the SFF archive.
 * It too is padded to an 8-byte boundary.
 *)
type sff_read_header = {
    hdr_len:int;
    name_len:int;
    nbases:int32;
    clip_qual_left:int;
    clip_qual_right:int;
    clip_adapter_left:int;
    clip_adapter_right:int;
    mutable name:string
  } 

(*
 * We have one read_data section per reading, following the read_header.
 * It is padded to an 8-byte boundary.
 *)
type sff_read_data = {
    flowgram:int list; (* x 100.0 *)
    flow_index:int list; (* relative to last *)
    bases:string;
    quality:int list
  }

type sffinfo = {
    sequence:string;
    signals:float list;
    quals:int list;
    indxs:int list
  } 

type ifht =
  | Ichan of in_channel
  | Iobj of sys_in_channel

type ofht =
  | Ochan of out_channel
  | Oobj of sys_out_channel

type read =
  | SFF454 of sffinfo
  | FASTQ of (string * string)
  | FASTQpe of read*read

let k = ref 0.15
let linsd1 = ref 0.1
let maxerror = ref 0.2
let insrate_Solexa = ref 0.0001
let delrate_Solexa = ref 0.0001

(* type of linker *)
type linker_type = NoLink | Titanium | Flx | Custom of string

let linker = ref NoLink

(* type paired = { linker_t:linker_type; *)
(*                 l_start:Int32; *)
(*                 l_length:Int32 *)
(*               }  *)

(* type r_info = { abbrev:string; *)
(*                 chrom:string; *)
(*                 from:int; *)
(*                 length:int; *)
(*                 str:char *)
(*               }  *)


type r_info = { length:int;
                str:char
              } 

type seqType = SES of string | PE454 of string | PESol of string * string
type headType = SEH of string | PEH of string * string
(*type headType = SEH of string | PEH of string * string*)
