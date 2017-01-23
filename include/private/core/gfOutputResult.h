/* gfOutputResult - Output controller for gfResult */

#ifndef _GFOUTPUTRESULT_H_
#define _GFOUTPUTRESULT_H_

/* free with freeGfOutputResult */
struct gfOutput *gfOutputResult(int minGood, boolean qIsProt, boolean tIsProt);

/* Free gfOutputResult */
void freeGfOutputResult(struct gfOutput **pp);

/* Free gfOutputResult but return the contained gfResult */
struct gfResult *unpackGfOutputResult(struct gfOutput **pp);

#endif /* GFOUTPUTRESULT_H */
