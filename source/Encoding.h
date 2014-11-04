/*

Pollux
Copyright (C) 2014  Eric Marinier

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

---

Contributions: Oliver Grant

*/

#ifndef ENCODING_H
#define	ENCODING_H

#ifdef	__cplusplus
extern "C" {
#endif
    
#include "Reads.h"

// Encodes the given string into the sequence_64b 64bit word array
void encode_sequence(struct read* rd, char* seq);

#ifdef	__cplusplus
}
#endif

#endif	/* ENCODING_H */
