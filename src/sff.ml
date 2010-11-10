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
open Types

(* decode_sff_common_header : string -> sff_common_header *)
let decode_sff_common_header str = 
  let strin = new string_in_channel str 0 in
  let magic = strin#read_int32 in
  let version = strin#read_string 4 in
  let index_offset = strin#read_int64 in
  let index_len = strin#read_int32 in
  let nreads = strin#read_int32 in
  let header_len = strin#read_int_size 2 in
  let key_len = strin#read_int_size 2 in
  let flow_len = strin#read_int_size 2 in
  let flowgram_format = strin#read_int_size 1 in
  { magic=magic;
    version=version;
    index_offset=index_offset;
    index_len=index_len;
    nreads=nreads;
    header_len=header_len;
    key_len=key_len;
    flow_len=flow_len;
    flowgram_format=flowgram_format;
    flow="";
    key=""
  } 

(* encode_sff_common_header : sff_common_header -> string *)
let encode_sff_common_header h =
  let buf = new buffer_out_channel (Buffer.create 1000) in
  begin
    buf#write_int32 h.magic;
    buf#write_string h.version;
    buf#write_int64 h.index_offset;
    buf#write_int32 h.index_len;
    buf#write_int32 h.nreads;
    buf#write_int16 h.header_len;
    buf#write_int16 h.key_len;
    buf#write_int16 h.flow_len;
    buf#write_byte  h.flowgram_format;
    buf#write_string h.flow;
    buf#write_string h.key;
    let e= 31 + h.flow_len + h.key_len in
    buf#write_string (String.sub (String.make 8 '\000') 0 ((((e+7) &! ~!7)-e)));
    buf#contents
  end

let read_sff_common_header cin =
  let chdr = cin#read_string 31 in
  let h = decode_sff_common_header chdr in
  begin
    h.flow <- cin#read_string h.flow_len;
    h.key <- cin#read_string h.key_len;
    let e = 31 + h.flow_len + h.key_len in
    ignore(cin#read_string ((((e + 7) &! ~!7)-e)));
    h
  end

let write_sff_common_header h cout =
  let str = encode_sff_common_header h in
  cout#write_string str
  
let decode_sff_read_header str =
  let cin = new string_in_channel str 0 in
  let hdr_len = cin#read_int_size 2 in
  let name_len = cin#read_int_size 2 in
  let nbases = cin#read_int32 in
  let clip_qual_left = cin#read_int_size 2 in
  let clip_qual_right = cin#read_int_size 2 in
  let clip_adapter_left = cin#read_int_size 2 in
  let clip_adapter_right = cin#read_int_size 2 in
  { hdr_len=hdr_len;
    name_len=name_len;
    nbases=nbases;
    clip_qual_left=clip_qual_left;
    clip_qual_right=clip_qual_right;
    clip_adapter_left=clip_adapter_left;
    clip_adapter_right=clip_adapter_right;
    name=""
  } 

let encode_sff_read_header h =
  let buf = new buffer_out_channel (Buffer.create 1000) in
  begin
    buf#write_int16 h.hdr_len;
    buf#write_int16 h.name_len;
    buf#write_int32 h.nbases;
    buf#write_int16 h.clip_qual_left;
    buf#write_int16 h.clip_qual_right;
    buf#write_int16 h.clip_adapter_left;
    buf#write_int16 h.clip_adapter_right;
    buf#write_string h.name;
    let e= 16 + h.name_len in
    buf#write_string (String.sub (String.make 8 '\000') 0 ((((e+7) &! ~!7)-e)));
    buf#contents
  end

let read_sff_read_header cin =
  let rhdr = cin#read_string 16 in
  let h = decode_sff_read_header rhdr in
  begin
    h.name <- cin#read_string h.name_len;
    let e = 16 + h.name_len in
    ignore(cin#read_string ((((e + 7) &! ~!7)-e)));
    h
  end

let write_nseqs n cout =
  begin
    cout#seek 20;
    cout#write_int32 (Int32.of_int n);
  end

let write_sff_read_header h cout =
  let str = encode_sff_read_header h in
  cout#write_string str

let read_sff_read_data cin nflows nbases =
  let flowgram = let arr = Array.make nflows 0 in
      begin
        for i = 0 to nflows - 1 do
          arr.(i) <- cin#read_int_size 2
        done;
        Array.to_list arr
      end in
  let flow_index = let arr = Array.make (Int32.to_int nbases) 0 in
      begin
        for i = 0 to ((Int32.to_int nbases) - 1) do
          arr.(i) <- let v = cin#read_int_size 1 in 
          try (v + arr.(i-1))
          with
          | Invalid_argument emsg ->
              if emsg = "index out of bounds" then v 
              else raise (Invalid_argument emsg)
        done;
        Array.to_list arr
      end in
  let bases = cin#read_string(Int32.to_int nbases) in
  let quality = let arr = Array.make (Int32.to_int nbases) 0 in
      begin
        for i = 0 to ((Int32.to_int nbases) - 1) do
          arr.(i) <- cin#read_int_size 1
        done;
        Array.to_list arr
      end in
  let h = 
    { flowgram=flowgram;
      flow_index=flow_index;
      bases=bases;
      quality=quality
    } in
  let e = 2*nflows + 3*(Int32.to_int nbases) in
  begin
    ignore(cin#read_string ((((e + 7) &! ~!7)-e)));
    h
  end

let write_sff_read_data h cout nflows nbases =
  begin
    List.iter (fun x -> cout#write_int16 x) h.flowgram;
    List.iter (fun x -> cout#write_byte x) h.flow_index;
    cout#write_string h.bases;
    List.iter (fun x -> cout#write_byte x) h.quality;
(* 8-byte padding *)
    let e = 2*nflows + 3*(Int32.to_int nbases) in
    cout#write_string (String.sub (String.make 8 '\000') 0 ((((e+7) &! ~!7)-e)))
  end

class insff fin =
  object(self)
    val cin = new sys_in_channel (open_in fin)
    val mutable common_header = {
      magic=sff_magic;
      version=sff_version;
      index_offset=Int64.zero;
      index_len=Int32.zero;
      nreads=Int32.zero;
      header_len=0;
      key_len=0;
      flow_len=0;
      flowgram_format=1;
      flow= "";
      key=""
    }
    method private read_sff_common_header = read_sff_common_header cin
    method private read_sff_read_header = read_sff_read_header cin
    method private read_sff_read_data nbases = read_sff_read_data cin common_header.flow_len nbases
    method readHeader = common_header
    method readRead =
      let rheader = self#read_sff_read_header in
      let rdata = self#read_sff_read_data rheader.nbases in
      (rheader,rdata)
    method close = cin#close
    initializer
      let cheader = self#read_sff_common_header in
      common_header <- cheader
  end

class outsff fout rounds =
  object(self)
    val cout = new sys_out_channel (open_out fout)
    val key = "TCAG"
    val mutable nseqs = 0
    val mutable ctrlseq = ""
    val flowOrder = "TACG"
    method write_sff_common_header h = write_sff_common_header h cout
    method write_sff_read_header h = write_sff_read_header h cout
    method write_sff_read_data h x y = 
      begin
        nseqs <- nseqs + 1;
        write_sff_read_data h cout x y
      end
    method writeRead_h (h,read) =
      match read with
      | SFF454 readRec ->
          let lseq = (String.length readRec.sequence) in 
          let r_header = {
            hdr_len = (let e= 16 + (String.length h) in ((e + 7) &! ~!7));
            name_len = (String.length h);
            nbases = Int32.of_int lseq;
            clip_qual_left = (String.length key)+1; 
            clip_qual_right = lseq;
            clip_adapter_left = 0;
            clip_adapter_right = 0;
            name = h
          } in
          let r_data = {
            flowgram = (List.map (fun x -> int_of_float (x *. 100.)) readRec.signals);
            flow_index = readRec.indxs;
            bases = readRec.sequence;
            quality = readRec.quals
          } in
          begin
            self#write_sff_read_header r_header;
            self#write_sff_read_data r_data (4*rounds) (Int32.of_int lseq)
          end
      | _ -> failwith "writeRead(454)(sff.ml) -- Invalid type of header or sequence" 
    method close = 
      begin
        write_nseqs nseqs cout;
        cout#close
      end
    initializer
      begin
        ctrlseq <- (
          let rec aux acc i = 
            if i >= rounds then acc 
            else aux (acc ^ flowOrder) (i+1)
          in aux "" 0
            );
        let c_header = {
          magic=sff_magic;
          version=sff_version;
          index_offset=Int64.zero;
          index_len=Int32.zero;
          nreads=Int32.of_int nseqs;
          header_len=(let e = 31 + (String.length key) + (rounds * 4) in ((e + 7) &! ~!7));
          key_len=(String.length key);
          flow_len=(rounds * 4);
          flowgram_format=1;
          flow= ctrlseq;
          key=key
        } in write_sff_common_header c_header cout
      end
  end
