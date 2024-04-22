from transformers import PreTrainedTokenizer, AddedToken
import collections
import os
from typing import List, Optional, Dict


""" AISTokenizer for Hugging Face Transformers.

"""
class AISTokenizer(PreTrainedTokenizer):
    def __init__(self, vocab: str, model_max_length: int, **kwargs):
        """Character tokens for Hugging Face transformers.

        Args:
            vocab str: Filename of a file containing with desired tokens
            on each newline
                    "<|endoftext|>": 0

                an id (starting at 1) will be assigned to each token.

            model_max_length (int): Model maximum sequence length.
        """
        self.vocab = []
        with open(vocab, 'r') as file:
            for line in file:
                self.vocab.append(line.strip())
        self.model_max_length = model_max_length
        
        bos_token = AddedToken("<|endoftext|>", lstrip=False, rstrip=False)
        eos_token = AddedToken("<|endoftext|>", lstrip=False, rstrip=False)
        pad_token = AddedToken("<|pad|>", lstrip=False, rstrip=False)
        unk_token = AddedToken("<|endoftext|>", lstrip=False, rstrip=False)
        self._vocab_str_to_int = {
            "<|endoftext|>": 0,"<|pad|>":1,
            **{ch: i + 2 for i, ch in enumerate(self.vocab)},
        }
        self._vocab_int_to_str = {v: k + " " for k, v in self._vocab_str_to_int.items()}
        self.ids_to_tokens = collections.OrderedDict([(ids, tok) for tok, ids in self._vocab_str_to_int.items()])

        super().__init__(
            unk_token=unk_token,
            bos_token=bos_token,
            eos_token=eos_token,
            pad_token=pad_token,
            add_prefix_space=False,
            model_max_length=model_max_length,
            **kwargs,
        )


    @property
    def vocab_size(self) -> int:
        return len(self._vocab_str_to_int)

    def _tokenize(self, text: str) -> List[str]:
        return text.split()

    def _convert_token_to_id(self, token: str) -> int:
        return self._vocab_str_to_int.get(token, self._vocab_str_to_int[self.unk_token])

    def _convert_id_to_token(self, index: int) -> str:
        return self._vocab_int_to_str[index]

    def convert_tokens_to_string(self, tokens):
        return "".join(tokens)

    def create_token_type_ids_from_sequences(
            self, token_ids_0: List[int], token_ids_1: Optional[List[int]] = None
        ) -> List[int]:
        bos_token_id = []
        eos_token_id = []

        output = [0] * len(bos_token_id + token_ids_0 + eos_token_id)

        if token_ids_1 is not None:
            output += [1] * len(bos_token_id + token_ids_1 + eos_token_id)

        return output
    def build_inputs_with_special_tokens(
        self, token_ids_0: List[int], token_ids_1: Optional[List[int]] = None
    ) -> List[int]:
        if True:
            bos_token_ids = [self.bos_token_id]
        else:
            bos_token_ids = []

        output = bos_token_ids + token_ids_0

        if token_ids_1 is None:
            return output

        return output + bos_token_ids + token_ids_1


    def get_special_tokens_mask(
        self,
        token_ids_0: List[int],
        token_ids_1: Optional[List[int]] = None,
        already_has_special_tokens: bool = False,
    ) -> List[int]:
        if already_has_special_tokens:
            if token_ids_1 is not None:
                raise ValueError(
                    "You should not supply a second sequence if the provided sequence of "
                    "ids is already formated with special tokens for the model."
                )
            return list(map(lambda x: 1 if x in [self.bos_token_id, self.eos_token_id] else 0, token_ids_0))

        if token_ids_1 is not None:
            return [1] + ([0] * len(token_ids_0)) + [1] + ([0] * len(token_ids_1)) + [1]
        return [1] + ([0] * len(token_ids_0)) + [1]
    def get_vocab(self) -> Dict[str, int]:
        return (self._vocab_str_to_int)
    def save_vocabulary(self, vocab_path,filename_prefix: Optional[str] = None):
        """
        Save the sentencepiece vocabulary (copy original file) and special tokens file to a directory.
        Args:
            vocab_path (:obj:`str`):
                The directory in which to save the vocabulary.
        Returns:
            :obj:`Tuple(str)`: Paths to the files saved.
        """
        index = 0
        if os.path.isdir(vocab_path):
            vocab_file = os.path.join(vocab_path, "vocab_file.txt")
        else:
            vocab_file = vocab_path
        with open(vocab_file, "w", encoding="utf-8") as writer:
            for token, token_index in sorted(self._vocab_str_to_int.items(), key=lambda kv: kv[1]):
                if index != token_index:
                    print(
                        "Saving vocabulary to {}: vocabulary indices are not consecutive."
                        " Please check that the vocabulary is not corrupted!".format(vocab_file)
                    )
                    index = token_index
                writer.write(token + "\n")
                index += 1
        return (vocab_file,)
