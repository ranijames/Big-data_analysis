import pandas as pd

HGNC_DB = '/home/mpschr/Data/hgnc/hgnc_symbols.txt'

__author__ = 'mpschr'


class Symbolizer(object):
    def __init__(self, hgnc_db_file=HGNC_DB,
                 keep_not_in_wanted=False,
                 update_unlisted=False
                 ):

        self._keep_not_in_wanted = keep_not_in_wanted
        self._update_unlisted = update_unlisted

        self.symbols = (
            pd.read_table(hgnc_db_file)
                .rename(columns=lambda x: x.replace(' ', '_').lower().replace('_symbol', ''))
        )

    def update_symbols(self, symbols, wanted_symbols: pd.Series=None):

        if type(symbols) == str:
            return self._update_symbol(symbols)

        # iterate over the symbol list
        try:
            iterator = iter(symbols)
        except TypeError:
            raise RuntimeError('The passed symbols must be iterable ()')

        if wanted_symbols is None:
            wanted_symbols = self.symbols.approved

        result = [self._update_symbol(x, wanted_symbols) for x in symbols]

        return result


    def _update_symbol(self, sym, wanted_symbols: pd.Series):

        sym = sym.strip()
        if sym in wanted_symbols.values:
            return sym

        approved = sym in self.symbols.approved.values
        previous = sym in self.symbols.previous.values
        alias = sym in self.symbols.alias.values

        foundas = [x for x, y in zip(['approved', 'previous', 'alias'], [approved, previous, alias]) if y]

        if approved and previous:
            print("Warning: {} is both previous and approved HGNC symbol".format(sym))

        if previous:
            updates = self.symbols[['approved', 'previous']].query('previous == "{}"'.format(sym)).drop_duplicates()
            if updates.shape[0] != 1:
                print(updates)
                raise RuntimeWarning("UPDATING PROBLEM for previous: " + sym)
            else:
                update = updates.iat[0, 0]
                if update not in wanted_symbols.values:
                    print(sym + ' failed - found in HGNC as: ', ','.join((foundas)))
                    if self._keep_not_in_wanted:
                        return sym
                    else:
                        return 'NOT_IN_WANTED;PREVIOUS'
                return update

        elif alias:
            updates = self.symbols[['approved', 'alias']].query('alias == "{}"'.format(sym)).drop_duplicates()
            update = None
            if updates.shape[0] == 1:
                update = updates.iat[0, 0]
            elif updates.shape[0] > 1:
                updates2 = (
                    self.symbols[['approved', 'previous', 'alias']]
                        .query('alias == "{}"'.format(sym)).drop_duplicates()
                        .fillna('')
                        .query('previous == ""')
                )
                print(updates2)
                if updates2.shape[0] == 1:
                    update = updates2.iat[0, 0]
                    print(
                        "Discarded {} for {}".format([x for x in updates.approved.drop_duplicates() if not x == update],
                                                     sym))

            if update is None:
                raise RuntimeWarning("UPDATING PROBLEM for alias: " + sym)
            else:
                update = updates.iat[0, 0]
                if update not in wanted_symbols.values:
                    print(sym + ' failed - found in HGNC as: ', ','.join((foundas)))
                    return 'NOT_IN_WANTED;ALIAS'
                return update
        else:
            print(sym + ' failed - found in HGNC as: ', ','.join((foundas)))
            if approved:
                if self._keep_not_in_wanted:
                    return sym
                else:
                    return 'NOT_IN_WANTED;APPROVED'
            else:
                if self._keep_not_in_wanted:
                    return sym
                else:
                    return 'NOT_IN_WANTED'
