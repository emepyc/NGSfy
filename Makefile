# Copyright 2010 Miguel Pignatelli. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

TARGET=ngsfy
all:
	cd src/myString && make && cd .. && make && cp $(TARGET) .. && cd ..

clean:
	cd src/myString && make clean && cd .. && make clean && cd ..
