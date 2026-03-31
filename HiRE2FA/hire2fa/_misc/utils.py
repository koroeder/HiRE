# //////////////////////////////////////////////////////////////////////////////
class Utils:
    # ------------------------------------------------------------------------------
    @staticmethod
    def argsort(arr: tuple | list) -> tuple[int]:
        """Helper function to get the indices that would sort a tuple/list"""
        if not len(arr): return ()
        idxs,_ = zip(*sorted(
            enumerate(arr), key = lambda t: int(t[1])
        ))
        return idxs


# //////////////////////////////////////////////////////////////////////////////
