#pragma once
#include "FullSystem/HessianBlocks.h"

namespace IOWrap
{
class Output3DWrapper
{
public:
    Output3DWrapper() {}
    virtual ~Output3DWrapper() {}

    /* call on finish */
    virtual void join() {}

    /* call on reset */
    virtual void reset() {}

    /* Usage:
     * Called once for each new frame, before it is tracked (i.e., it doesn't have a pose yet).
     *
     * Calling:
     * Always called, no overhead if not used.
     */
    virtual void pushLiveFrame(FrameHessian* image) {}
};
}