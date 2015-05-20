#ifndef ZEROGRADIENTPATCH_H
#define ZEROGRADIENTPATCH_H


#ifndef FIXEDGRADIENTPATCH_H
#define FIXEDGRADIENTPATCH_H

#include "PatchField.h"

namespace SEM { namespace field {


//TODO - PatchField need to be modified. Some implementaion
//       have to go to BasicFiledPatch.
//       It's neccessery due to fact, that ZeroGradient patch
//       do not need to store any values, it's some kind of
//       dummy patch

}//field
}//SEM


#endif // FIXEDGRADIENTPATCH_H


#endif // ZEROGRADIENTPATCH_H
