'''
This tool combs through the code and keeps a tally of how many indentations
each line of code should have (and automatically corrects it). Not too much to
say about it besides the logic can get a little twisted. Each line is checked
twice; once before adjusting whitespace, and once after (to affect the next line).

Note: this is the only component to print errors to standard output. Two errors
can occur, and both involve the `self.levelNumber` variable. This tracks how many
blocks in the code should be. For example, being inside a subroutine will increment
self.levelNumber by 1, as will being inside an if statement.

At the end of the program, self.levelNumber should be 0. If not, this means there
is a code block that is not terminated properly, or perhaps a bug in the linter's
ability to detect them.

The other error is thrown if self.levelNumber ever drops below 0, which should not
happen. If this error arises, check the line number where it happens.
'''
import re

class WhitespaceChecker(object):

    def __init__(self, fname):
        self.fname = fname
        with open(self.fname, 'r') as f:
            lines = f.readlines()
        self.lines = lines
        self.fixedlines = []
        self.ignoreLines = 0
        self.directive = False
        self.levelNumber = 0
        self.continuedIf = False
        self.checkWhiteSpace()
        if self.levelNumber != 0:
            print("Whitespace error in {}: finished with levelNumber {}".format(self.fname, self.levelNumber))
        with open(self.fname, 'w') as f:
            for line in self.fixedlines:
                if not line.endswith('\n'):
                    line += '\n'
                f.write(line)

    def checkWhiteSpace(self):
        stack = []
        continuationLine = False
        subcall = 0
        errr = False
        for i, line in enumerate(self.lines):
            if self.ignoreLines == 0 and self.levelNumber < 0 and not errr:
                print("Whitepsace error in {} at line {}: levelNumber < 0".format(self.fname, i+1))
                errr = True
            # ignore blank lines
            if re.match(r'^\s*\n', line, re.IGNORECASE):
                self.fixedlines.append("\n")
                continue
            if self.skippable(line):
                continue
            elif continuationLine:
                newline = ""
                if subcall > 0:
                    for x in range(subcall):
                        newline += " "
                    newline += line.strip()
                else:
                    newline = line
                # end of continuation lines
                if not re.match(r'^[^\n!]*&', line, re.IGNORECASE) and \
                   not re.match(r'^\s*!', line, re.IGNORECASE):
                    continuationLine = False
                    subcall = 0
                self.fixedlines.append(newline)
                continue
            else:
                newline = ""
                # check that this line follows regular formatting
                parsed_line = re.match(r'^(\s*)(\d*)(\s*)(.*)\n', line, re.IGNORECASE)
                if not parsed_line:
                    self.fixedlines.append(line)
                    continue

                self.analyzeLine1(parsed_line.group(4))

                # if we have a labelled line, we must account for that in the whitespace
                if parsed_line.group(2) != "":
                    newline += parsed_line.group(2)
                    newline += " "

                # add whitespace based on the levelNumber, and account for line labelling
                num_spaces = max(0, 2*self.levelNumber - len(newline))
                for x in range(num_spaces):
                    newline += " "
                if re.match(r'\bcontains\b', parsed_line.group(4), re.IGNORECASE):
                    newline = newline[:-2]
                newline += parsed_line.group(4)
                self.fixedlines.append(newline)
                self.analyzeLine2(parsed_line.group(4))
                if re.match(r'^[^\n!]*&', line, re.IGNORECASE) and not self.continuedIf:
                    continuationLine = True
                    line2 = re.match(r'^([^\n!]*)\b(call)\b([^\n!\(]*)\(', line, re.IGNORECASE)
                    if line2:
                        subcall = len(line2.group(0))

    # checks if a line should be skipped by the linter
    def skippable(self, line):
        if self.directive:
            if re.match(r'^[ \t]*#endif', line, re.IGNORECASE):
                self.directive = False
            self.fixedlines.append(line)
            return True
        else:
            if re.match(r'^[ \t]*#if', line, re.IGNORECASE):
                self.fixedlines.append(line)
                self.directive = True
                return True
            else:
                lineSkip = re.match(r'^\s*!ignoreLint\((\d+)\)ff', line, re.IGNORECASE)
                if lineSkip:
                    self.ignoreLines = int(lineSkip.group(1))
                    self.fixedlines.append(line)
                    return True
                elif self.ignoreLines > 0:
                    self.ignoreLines -= 1
                    self.fixedlines.append(line)
                    return True
        return False

    # checks if we should revert whitespace on this line (eg. 'end if')
    def analyzeLine1(self, line):
        if re.match(r'^[^!]*<<<', line, re.IGNORECASE): # ignore lines that the linter has already flagged as awful
            return
        elif re.match(r'^\s*\b(end|case|else)\b', line, re.IGNORECASE):
            self.levelNumber -= 1

    # checks if we should increase whitespace on the NEXT line (eg. 'do j=1,n')
    def analyzeLine2(self, line):
        if re.match(r'^[^!]*<<<', line, re.IGNORECASE): # ignore lines that the linter has already flagged as awful
            return
        elif re.match(r'^[^\n!\'"]*\b(module(?! procedure )|interface|function|program|type |do|subroutine|case|select|else|forall|where)\b', line, re.IGNORECASE) \
        and not re.match(r'^[^\n!]*\bend\b', line, re.IGNORECASE):
            self.levelNumber += 1
            self.continuedIf = False
        elif re.match(r'^[^!\n\'"]*\bif\b[^!\n]*\bthen\b', line, re.IGNORECASE):
            self.levelNumber += 1
            self.continuedIf = False
        elif re.match(r'^[^!\n\'"]*(?<!end )\bif\b[^!\n]*&', line, re.IGNORECASE):
            self.continuedIf = True
        elif self.continuedIf and re.match(r'^[^!\n]*\bthen\b', line, re.IGNORECASE):
            self.continuedIf = False
            self.levelNumber += 1
